#include "IndexCreator.h"
#include "FileUtil.h"
#include "LocalUtil.h"
#include "ProdigalWrapper.h"
#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <utility>
#include "NcbiTaxonomy.cpp"
#include "common.h"
#include "geneInfoAdder.h"

extern const char *version;

geneInfoAdder::geneInfoAdder(const LocalParameters &par, TaxonomyWrapper *taxonomy) : IndexCreator(par, taxonomy)
{
  dbDir = par.filenames[0];
  if (par.taxonomyPath.empty())
  {
    taxonomyDir = dbDir + "/taxonomy/";
  }
  else
  {
    taxonomyDir = par.taxonomyPath + "/";
  }
  fnaListFileName = par.filenames[1];
  acc2taxidFileName = par.filenames[2];

  taxidListFileName = dbDir + "/taxID_list";
  geneListFileName = dbDir + "/geneID_list";
  taxonomyBinaryFileName = dbDir + "/taxonomyDB";
  versionFileName = dbDir + "/db.version";
  paramterFileName = dbDir + "/db.parameters";

  if (par.reducedAA == 1)
  {
    MARKER = 0Xffffffff;
    MARKER = ~MARKER;
  }
  else
  {
    MARKER = 16777215;
    MARKER = ~MARKER;
  }

  isUpdating = false;
  subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
}

geneInfoAdder::~geneInfoAdder()
{
  delete subMat;
}
void geneInfoAdder::createGeneIndex(const LocalParameters &par) {
    TargetKmerBuffer kmerBuffer(calculateBufferSize(par.ramUsage));
    indexReferenceSequences(kmerBuffer.bufferSize);
    
    if (!par.cdsInfo.empty()) {
        loadCdsInfo(par.cdsInfo);   
    }

    cout << "Made blocks for each thread" << endl;

    // Write taxonomy id list
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid: taxIdSet) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);

    // Process the splits until all are processed
    size_t batchNum = accessionBatches.size();
    bool * batchChecker = new bool[batchNum];
    fill_n(batchChecker, batchNum, false);
    size_t processedBatchCnt = 0;
    
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    cout << "Kmer buffer size: " << kmerBuffer.bufferSize << endl;
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(processedBatchCnt < batchNum) {
        memset(kmerBuffer.buffer, 0, kmerBuffer.bufferSize * sizeof(TargetKmer));
        fillTargetKmerBuffer(kmerBuffer, batchChecker, processedBatchCnt, par);

        // Sort the k-mers
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      IndexCreator::compareForDiffIdx);
        time_t sort = time(nullptr);
        cout << "Sort time: " << sort - start << endl;

        // Reduce redundancy
        auto * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        reduceRedundancy(kmerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges, par);
        time_t reduction = time(nullptr);
        cout << "Time spent for reducing redundancy: " << (double) (reduction - sort) << endl;
        if(processedBatchCnt == batchNum && numOfFlush == 0 && !isUpdating) {
            writeTargetFilesAndSplits(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        } else {
            writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, uniqKmerIdx, uniqKmerIdxRanges);
        }
        delete[] uniqKmerIdx;
    }
    delete[] batchChecker;
    writeDbParameters();
}

size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer &kmerBuffer,
                                          bool *checker,
                                          size_t &processedBatchCnt,
                                          const LocalParameters &par) {
    int hasOverflow = 0;
#pragma omp parallel default(none), shared(kmerBuffer, checker, processedBatchCnt, hasOverflow, par, cout)
    {
        ProbabilityMatrix probMatrix(*subMat);
        // ProdigalWrapper * prodigal = new ProdigalWrapper();
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t orfNum;
        vector<PredictedBlock> extendedORFs;
        priority_queue<uint64_t> standardList;
        priority_queue<uint64_t> currentList;
        size_t lengthOfTrainingSeq;
        char *reverseCompliment;
        vector<uint64_t> intergenicKmers;
        vector<int> aaSeq;
        vector<string> cds;
        vector<string> nonCds;
        bool trained = false;
		
#pragma omp for schedule(dynamic, 1)
        for (size_t batchIdx = 0; batchIdx < accessionBatches.size(); batchIdx ++) {
            if (!checker[batchIdx] && !hasOverflow) {
                checker[batchIdx] = true;
                intergenicKmers.clear();
                standardList = priority_queue<uint64_t>();

                // Estimate the number of k-mers to be extracted from current split
                size_t totalLength = 0;
                for (size_t p = 0; p < accessionBatches[batchIdx].lengths.size(); p++) {
                    totalLength += accessionBatches[batchIdx].lengths[p];
                }
                size_t estimatedKmerCnt = (totalLength + totalLength / 5) / 3;

                ProdigalWrapper * prodigal = new ProdigalWrapper();
                trained = false;

                // Process current split if buffer has enough space.
                posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
                if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize) {
                    KSeqWrapper* kseq = KSeqFactory(fastaPaths[accessionBatches[batchIdx].whichFasta].c_str());
                    #ifdef NDEBUG
                        #pragma omp critical
                        {
                            cout << fastaPaths[accessionBatches[batchIdx].whichFasta] << endl; 
                        }
                    #endif
                    size_t seqCnt = 0;
                    size_t idx = 0;
                    while (kseq->ReadEntry()) {
                        if (seqCnt == accessionBatches[batchIdx].orders[idx]) {
                            const KSeqWrapper::KSeqEntry & e = kseq->entry;
                            #ifdef NDEBUG
                                #pragma omp critical
                                {
                                    cout << "Processing " << e.name.s << "\t" 
                                         << e.sequence.l << "\t" 
                                         << taxonomy->getOriginalTaxID(accessionBatches[batchIdx].speciesID) << "\t" 
                                         << taxonomy->getOriginalTaxID(accessionBatches[batchIdx].taxIDs[idx]) << "\n";
                                }
                            #endif

                            // Mask low complexity regions
                            char *maskedSeq = nullptr;
                            if (par.maskMode) {
                                maskedSeq = new char[e.sequence.l + 1]; // TODO: reuse the buffer
                                SeqIterator::maskLowComplexityRegions((unsigned char *) e.sequence.s, (unsigned char *) maskedSeq, probMatrix, par.maskProb, subMat);
                                maskedSeq[e.sequence.l] = '\0';
                            } else {
                                maskedSeq = e.sequence.s;
                            }

                            orfNum = 0;
                            extendedORFs.clear();
                            int tempCheck = 0;                            
                            if (cdsInfoMap.find(string(e.name.s)) != cdsInfoMap.end()) {
                                // Get CDS and non-CDS
                                cds.clear();
                                nonCds.clear();
                                seqIterator.devideToCdsAndNonCds(maskedSeq,
                                                                 e.sequence.l,
                                                                 cdsInfoMap[string(e.name.s)],
                                                                 cds,
                                                                 nonCds);
    
                                for (size_t cdsCnt = 0; cdsCnt < cds.size(); cdsCnt ++) {
                                    seqIterator.translate(cds[cdsCnt], aaSeq);
                                    tempCheck = seqIterator.fillBufferWithKmerFromBlock(
                                                    cds[cdsCnt].c_str(),
                                                    kmerBuffer,
                                                    posToWrite,
                                                    accessionBatches[batchIdx].taxIDs[idx],
                                                    accessionBatches[batchIdx].speciesID,
                                                    aaSeq);
                                    if (tempCheck == -1) {
                                        cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                    }
                                }

                                for (size_t nonCdsCnt = 0; nonCdsCnt < nonCds.size(); nonCdsCnt ++) {
                                    seqIterator.translate(nonCds[nonCdsCnt], aaSeq);
                                    tempCheck = seqIterator.fillBufferWithKmerFromBlock(
                                                    nonCds[nonCdsCnt].c_str(),
                                                    kmerBuffer,
                                                    posToWrite,
                                                    accessionBatches[batchIdx].taxIDs[idx],
                                                    accessionBatches[batchIdx].speciesID,
                                                    aaSeq);
                                    if (tempCheck == -1) {
                                        cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                    }
                                }
                            } else {
                                // USE PRODIGAL
                                if (!trained) {
                                    KSeqWrapper* training_seq = KSeqFactory(fastaPaths[accessionBatches[batchIdx].trainingSeqFasta].c_str()); 
                                    size_t seqCnt = 0;
                                    while (training_seq->ReadEntry()) {
                                        if (seqCnt == accessionBatches[batchIdx].trainingSeqIdx) {
                                            break;
                                        }
                                        seqCnt++;
                                    }
                                    // #pragma omp critical
                                    // {
                                    //     cout << "Training " << accessionBatches[batchIdx].speciesID << "\t" << training_seq->entry.name.s << "\t" << training_seq->entry.sequence.l << endl;
                                    // }
                                    lengthOfTrainingSeq = training_seq->entry.sequence.l;
                                    prodigal->is_meta = 0;
                                    if (lengthOfTrainingSeq < 100'000) {
                                        prodigal->is_meta = 1;
                                        prodigal->trainMeta((unsigned char *) training_seq->entry.sequence.s, training_seq->entry.sequence.l);
                                    } else {
                                        prodigal->trainASpecies((unsigned char *) training_seq->entry.sequence.s, training_seq->entry.sequence.l);
                                    }

                                    // Generate intergenic 23-mer list. It is used to determine extension direction of intergenic sequences.
                                    prodigal->getPredictedGenes((unsigned char *) training_seq->entry.sequence.s, training_seq->entry.sequence.l);
                                    seqIterator.generateIntergenicKmerList(prodigal->genes, prodigal->nodes,
                                                                           prodigal->getNumberOfPredictedGenes(),
                                                                           intergenicKmers, training_seq->entry.sequence.s);

                                    // Get min k-mer hash list for determining strandness
                                    seqIterator.getMinHashList(standardList, training_seq->entry.sequence.s);
                                    delete training_seq;
                                    trained = true;
                                }

                                currentList = priority_queue<uint64_t>();
                                seqIterator.getMinHashList(currentList, e.sequence.s);

                                if (seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq, e.sequence.l)) {
                                    // Get extended ORFs
                                    prodigal->getPredictedGenes((unsigned char *) e.sequence.s, e.sequence.l);
                                    prodigal->removeCompletelyOverlappingGenes();
                                    prodigal->getExtendedORFs(prodigal->finalGenes, prodigal->nodes, extendedORFs,
                                                                 prodigal->fng, e.sequence.l,
                                                            orfNum, intergenicKmers, e.sequence.s);

                                    // Get k-mers from extended ORFs
                                    for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                        aaSeq.clear();
                                        seqIterator.translateBlock(maskedSeq, extendedORFs[orfCnt], aaSeq, e.sequence.l);
                                        tempCheck = seqIterator.fillBufferWithKmerFromBlock(
                                                extendedORFs[orfCnt],
                                                maskedSeq,
                                                kmerBuffer,
                                                posToWrite,
                                                accessionBatches[batchIdx].taxIDs[idx],
                                                accessionBatches[batchIdx].speciesID,
                                                aaSeq);
                                        if (tempCheck == -1) {
                                            cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                        }
                                    }
                                } else { // Reverse complement
                                    reverseCompliment = seqIterator.reverseCompliment(e.sequence.s, e.sequence.l);

                                    // Get extended ORFs
                                    prodigal->getPredictedGenes((unsigned char *) reverseCompliment, e.sequence.l);
                                    prodigal->removeCompletelyOverlappingGenes();
                                    prodigal->getExtendedORFs(prodigal->finalGenes, prodigal->nodes, extendedORFs,
                                                                     prodigal->fng, e.sequence.l,
                                                                orfNum, intergenicKmers, reverseCompliment);

                                    // Get reverse masked sequence
                                    if (par.maskMode) {
                                        delete[] maskedSeq;
                                        maskedSeq = new char[e.sequence.l + 1];
                                        SeqIterator::maskLowComplexityRegions((unsigned char *) reverseCompliment, (unsigned char *) maskedSeq, probMatrix, par.maskProb, subMat);
                                        maskedSeq[e.sequence.l] = '\0';
                                    } else {
                                        maskedSeq = reverseCompliment;
                                    }

                                    for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                        aaSeq.clear();
                                        seqIterator.translateBlock(maskedSeq, extendedORFs[orfCnt], aaSeq, e.sequence.l);
                                        tempCheck = seqIterator.fillBufferWithKmerFromBlock(
                                                extendedORFs[orfCnt],
                                                maskedSeq,
                                                kmerBuffer,
                                                posToWrite,
                                                accessionBatches[batchIdx].taxIDs[idx],
                                                accessionBatches[batchIdx].speciesID,
                                                aaSeq);
                                        if (tempCheck == -1) {
                                            cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                        }
                                    }
                                    free(reverseCompliment);  
                                }                            
                            }
                            idx++;
                            if (par.maskMode) {
                                delete[] maskedSeq;
                            }
                            if (idx == accessionBatches[batchIdx].lengths.size()) {
                                break;
                            }
                        }
                        seqCnt++;
                    }
                    delete kseq;
                    __sync_fetch_and_add(&processedBatchCnt, 1);
                    #pragma omp critical
                    {
                        cout << processedBatchCnt << " batches processed out of " << accessionBatches.size() << endl;
                    }
                } else {
                    checker[batchIdx] = false;
                    __sync_fetch_and_add(&hasOverflow, 1);
                    __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
                }
                delete prodigal;   
            }
        }
        // delete prodigal;
    }

    cout << "Before return: " << kmerBuffer.startIndexOfReserve << endl;
    return 0;
}

