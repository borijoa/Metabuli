#include "ProdigalWrapper.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "omp.h"
#include <vector>
#include <iostream>
#include <string>
#include <dirent.h>
#include <KSeqWrapper.h>
#include <sys/stat.h>
#include <sys/types.h>

size_t getSequenceLength(fptr fp);
bool hasReadPermission(const std::string &filePath);

int prodigal(const std::string &dbDir, const LocalParameters &par);

int prodigal(int argc, const char **argv, const Command &command)
{
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string &dbDir = par.filenames[0];

    if (!FileUtil::directoryExists(dbDir.c_str()))
    {
        Debug(Debug::INFO) << "DB directory" << dbDir << " is NOT exists.\n";
        return 0;
    }

    return prodigal(dbDir, par);
}

void ensureDirectoryExists(const std::string &path)
{
    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR))
    {
        Debug(Debug::INFO) << "Creating directory: " << path << "\n";
        if (mkdir(path.c_str(), 0777) != 0)
        {
            Debug(Debug::ERROR) << "Failed to create directory: " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
};

int prodigal(const std::string &dbDir, const LocalParameters &par)
{

    // library and cds path options

    const string libraryPath = dbDir + "library/";
    const string cdsPath = dbDir + "cds/";
    const string aaPath = dbDir + "aa/";
    ensureDirectoryExists(cdsPath);
    ensureDirectoryExists(aaPath);

    DIR *dir = opendir(libraryPath.c_str());
    if (dir == nullptr)
    {
        Debug(Debug::ERROR) << "Could not open " << libraryPath << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    
        // ...existing code...

        std::vector<std::pair<std::string, off_t>> fnaFilesWithSize;

        struct dirent *entry;
        while ((entry = readdir(dir)) != nullptr)
        {
            // Skip "." and ".."
            if (std::string(entry->d_name) == "." || std::string(entry->d_name) == "..")
            {
                continue;
            }

            // Check for .fna extension
            std::string fileName = entry->d_name;
            if (fileName.size() > 4 && fileName.substr(fileName.size() - 4) == ".fna")
            {
                std::string fullPath = libraryPath + fileName;

                // Get file size
                struct stat fileStat;
                if (stat(fullPath.c_str(), &fileStat) == 0)
                {
                    fnaFilesWithSize.emplace_back(fullPath, fileStat.st_size);
                }
                else
                {
                    Debug(Debug::ERROR) << "Could not get size for file: " << fullPath << "\n";
                }
            }
        }
        closedir(dir);

        // Sort files by size
        std::sort(fnaFilesWithSize.begin(), fnaFilesWithSize.end(),
                  [](const std::pair<std::string, off_t> &a, const std::pair<std::string, off_t> &b)
                  {
                      return a.second < b.second;
                  });

        // Extract sorted file paths
        std::vector<std::string> fnaFiles;
        
        for (const auto &file : fnaFilesWithSize)
        {
            fnaFiles.push_back(file.first);
        }
        std::vector<bool> fileProcessed(fnaFiles.size(), false);
        omp_set_num_threads(par.threads);

#pragma omp parallel default(none) shared(stdout, stderr, fnaFiles, cerr, cout, libraryPath, cdsPath, aaPath, par, fileProcessed)
    {
        

#pragma omp for schedule(static, 1)
        for (size_t a = 0; a < fnaFiles.size(); ++a)
        {
            #pragma omp critical
            {
                std::cout << "Thread " << omp_get_thread_num() << " processing: " << fnaFiles[a] << std::endl;
            }
        
            std::string cdsFileName = cdsPath + fnaFiles[a].substr(libraryPath.size(), fnaFiles[a].size() - libraryPath.size() - 4) + "_cds";
            std::string aaFileName = aaPath + fnaFiles[a].substr(libraryPath.size(), fnaFiles[a].size() - libraryPath.size() - 4) + "_aa";

            // FILE *inputFile = fopen(fnaFiles[a].c_str(), "r");
            FILE *cdsFile = fopen(cdsFileName.c_str(), "w");
            FILE *aaFile = fopen(aaFileName.c_str(), "w");
            FILE *inputFile = fopen(fnaFiles[a].c_str(), "r");

            int rv, slen, nn, ng, i, ipath, *gc_frame, do_training, output, max_phase;
            int closed, do_mask, nmask, force_nonsd, user_tt, is_meta, num_seq, quiet;
            int piped, max_slen, fnum;
            double max_score, gc, low, high;
            unsigned char *seq, *rseq, *useq;
            char *train_file, *start_file, *trans_file, *nuc_file;
            char *input_file, *output_file, input_copy[MAX_LINE];
            char cur_header[MAX_LINE], new_header[MAX_LINE], short_header[MAX_LINE];
            FILE *output_ptr, *start_ptr, *trans_ptr, *nuc_ptr;
            fptr input_ptr = inputFile;
            struct stat fbuf;
            pid_t pid;
            struct _node *nodes;
            struct _gene *genes;
            struct _training tinf;
            struct _metagenomic_bin meta[NUM_META];
            mask mlist[MAX_MASKS];

            /* Allocate memory and initialize variables */
            seq = (unsigned char *)malloc(MAX_SEQ / 4 * sizeof(unsigned char));
            rseq = (unsigned char *)malloc(MAX_SEQ / 4 * sizeof(unsigned char));
            useq = (unsigned char *)malloc(MAX_SEQ / 8 * sizeof(unsigned char));
            nodes = (struct _node *)malloc(STT_NOD * sizeof(struct _node));
            genes = (struct _gene *)malloc(MAX_GENES * sizeof(struct _gene));
            if (seq == NULL || rseq == NULL || nodes == NULL || genes == NULL)
            {
                fprintf(stderr, "\nError: Malloc failed on sequence/orfs\n\n");
                exit(1);
            }
            memset(seq, 0, MAX_SEQ / 4 * sizeof(unsigned char));
            memset(rseq, 0, MAX_SEQ / 4 * sizeof(unsigned char));
            memset(useq, 0, MAX_SEQ / 8 * sizeof(unsigned char));
            memset(nodes, 0, STT_NOD * sizeof(struct _node));
            memset(genes, 0, MAX_GENES * sizeof(struct _gene));
            memset(&tinf, 0, sizeof(struct _training));

            for (int i = 0; i < NUM_META; i++)
            {
                memset(&meta[i], 0, sizeof(struct _metagenomic_bin));
                strcpy(meta[i].desc, "None");
                meta[i].tinf = (struct _training *)malloc(sizeof(struct _training));
                if (meta[i].tinf == NULL)
                {
                    fprintf(stderr, "\nError: Malloc failed on training structure.\n\n");
                    exit(1);
                }
                memset(meta[i].tinf, 0, sizeof(struct _training));
            }
            nn = 0;
            slen = 0;
            ipath = 0;
            ng = 0;
            nmask = 0;
            user_tt = 0;
            is_meta = 0;
            num_seq = 0;
            quiet = 1;
            max_phase = 0;
            max_score = -100.0;
            train_file = NULL;
            do_training = 0;
            start_file = NULL;
            trans_file = NULL;
            nuc_file = NULL;
            start_ptr = stdout;
            trans_ptr = aaFile;
            nuc_ptr = stdout;
            input_file = NULL;
            output_file = NULL;
            piped = 0;
            output_ptr = cdsFile;
            max_slen = 0;
            output = 3;
            closed = 1;
            do_mask = 0;
            force_nonsd = 0;

            /* Filename for input copy if needed */
            pid = getpid();
            sprintf(input_copy, "tmp.prodigal.stdin.%d", pid);

            /***************************************************************************
              Set the start score weight.  Changing this number can dramatically
              affect the performance of the program.  Some genomes want it high (6+),
              and some prefer it low (2.5-3).  Attempts were made to determine this
              weight dynamically, but none were successful.  Therefore, we just
              manually set the weight to an average value that seems to work decently
              for 99% of genomes.  This problem may be revisited in future versions.
            ***************************************************************************/
            tinf.st_wt = 4.35;
            tinf.trans_table = 11;

            size_t sLen = getSequenceLength(input_ptr);

            // choose meta or single with sequence length
            if (sLen < 100'000)
            {
                is_meta = 1;
                cout << "META" << endl;
            }

            rewind(input_ptr);

            if (is_meta == 0 && (do_training == 1 || (do_training == 0 && train_file ==
                                                                              NULL)))
            {
                if (quiet == 0)
                {
                    fprintf(stderr, "Request:  Single Genome, Phase:  Training\n");
                    fprintf(stderr, "Reading in the sequence(s) to train...");
                }
                slen = read_seq_training(input_ptr, seq, useq, &(tinf.gc), do_mask, mlist,
                                         &nmask);
                if (slen == 0)
                {
                    fprintf(stderr, "\n\nSequence read failed (file must be Fasta, ");
                    fprintf(stderr, "Genbank, or EMBL format).\n\n");
                    exit(9);
                }
                if (slen < MIN_SINGLE_GENOME)
                {
                    fprintf(stderr, "\n\nError:  Sequence must be %d", MIN_SINGLE_GENOME);
                    fprintf(stderr, " characters (only %d read).\n(Consider", slen);
                    fprintf(stderr, " running with the -p meta option or finding");
                    fprintf(stderr, " more contigs from the same genome.)\n\n");
                    exit(10);
                }
                if (slen < IDEAL_SINGLE_GENOME)
                {
                    fprintf(stderr, "\n\nWarning:  ideally Prodigal should be given at");
                    fprintf(stderr, " least %d bases for ", IDEAL_SINGLE_GENOME);
                    fprintf(stderr, "training.\nYou may get better results with the ");
                    fprintf(stderr, "-p meta option.\n\n");
                }
                rcom_seq(seq, rseq, useq, slen);
                if (quiet == 0)
                {
                    fprintf(stderr, "%d bp seq created, %.2f pct GC\n", slen, tinf.gc * 100.0);
                }

                /***********************************************************************
                  Find all the potential starts and stops, sort them, and create a
                  comprehensive list of nodes for dynamic programming.
                ***********************************************************************/
                if (quiet == 0)
                {
                    fprintf(stderr, "Locating all potential starts and stops...");
                }
                if (slen > max_slen && slen > STT_NOD * 8)
                {
                    nodes = (struct _node *)realloc(nodes, (int)(slen / 8) * sizeof(struct _node));
                    if (nodes == NULL)
                    {
                        fprintf(stderr, "Realloc failed on nodes\n\n");
                        exit(11);
                    }
                    max_slen = slen;
                }

                nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
                qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
                if (quiet == 0)
                {
                    fprintf(stderr, "%d nodes\n", nn);
                }

                /***********************************************************************
                  Scan all the ORFS looking for a potential GC bias in a particular
                  codon position.  This information will be used to acquire a good
                  initial set of genes.
                ***********************************************************************/
                if (quiet == 0)
                {
                    fprintf(stderr, "Looking for GC bias in different frames...");
                }
                gc_frame = calc_most_gc_frame(seq, slen);
                if (gc_frame == NULL)
                {
                    fprintf(stderr, "Malloc failed on gc frame plot\n\n");
                    exit(11);
                }
                record_gc_bias(gc_frame, nodes, nn, &tinf);
                if (quiet == 0)
                {
                    fprintf(stderr, "frame bias scores: %.2f %.2f %.2f\n", tinf.bias[0],
                            tinf.bias[1], tinf.bias[2]);
                }
                free(gc_frame);

                /***********************************************************************
                  Do an initial dynamic programming routine with just the GC frame
                  bias used as a scoring function.  This will get an initial set of
                  genes to train on.
                ***********************************************************************/
                if (quiet == 0)
                {
                    fprintf(stderr, "Building initial set of genes to train from...");
                }
                record_overlapping_starts(nodes, nn, &tinf, 0);
                ipath = dprog(nodes, nn, &tinf, 0);
                if (quiet == 0)
                {
                    fprintf(stderr, "done!\n");
                }

                /***********************************************************************
                  Gather dicodon statistics for the training set.  Score the entire set
                  of nodes.
                ***********************************************************************/
                if (quiet == 0)
                {
                    fprintf(stderr, "Creating coding model and scoring nodes...");
                }
                calc_dicodon_gene(&tinf, seq, rseq, slen, nodes, ipath);
                raw_coding_score(seq, rseq, slen, nodes, nn, &tinf);
                if (quiet == 0)
                {
                    fprintf(stderr, "done!\n");
                }

                /***********************************************************************
                  Determine if this organism uses Shine-Dalgarno or not and score the
                  nodes appropriately.
                ***********************************************************************/
                if (quiet == 0)
                {
                    fprintf(stderr, "Examining upstream regions and training starts...");
                }
                rbs_score(seq, rseq, slen, nodes, nn, &tinf);
                train_starts_sd(seq, rseq, slen, nodes, nn, &tinf);
                determine_sd_usage(&tinf);
                if (force_nonsd == 1)
                    tinf.uses_sd = 0;
                if (tinf.uses_sd == 0)
                    train_starts_nonsd(seq, rseq, slen, nodes, nn, &tinf);
                if (quiet == 0)
                {
                    fprintf(stderr, "done!\n");
                }

                /* If training specified, write the training file and exit. */
                if (do_training == 1)
                {
                    if (quiet == 0)
                    {
                        fprintf(stderr, "Writing data to training file %s...", train_file);
                    }
                    rv = write_training_file(train_file, &tinf);
                    if (rv != 0)
                    {
                        fprintf(stderr, "\nError: could not write training file!\n");
                        exit(12);
                    }
                    else
                    {
                        if (quiet == 0)
                            fprintf(stderr, "done!\n");
                        exit(0);
                    }
                }

                /* Rewind input file */
                if (quiet == 0)
                    fprintf(stderr, "-------------------------------------\n");
                rewind(input_ptr);
                if (INPUT_SEEK(input_ptr, 0, SEEK_SET) == -1)
                {
                    fprintf(stderr, "\nError: could not rewind input file.\n");
                    exit(13);
                }

                /* Reset all the sequence/dynamic programming variables */
                memset(seq, 0, (slen / 4 + 1) * sizeof(unsigned char));
                memset(rseq, 0, (slen / 4 + 1) * sizeof(unsigned char));
                memset(useq, 0, (slen / 8 + 1) * sizeof(unsigned char));
                memset(nodes, 0, nn * sizeof(struct _node));
                nn = 0;
                slen = 0;
                ipath = 0;
                nmask = 0;
            }

            else if (is_meta == 1)
            {
                if (quiet == 0)
                {
                    fprintf(stderr, "Request:  Metagenomic, Phase:  Training\n");
                    fprintf(stderr, "Initializing training files...");
                }
                initialize_metagenomic_bins(meta);
                if (quiet == 0)
                {
                    fprintf(stderr, "done!\n");
                    fprintf(stderr, "-------------------------------------\n");
                }
            }

            /* Print out header for gene finding phase */
            if (quiet == 0)
            {
                if (is_meta == 1)
                    fprintf(stderr, "Request:  Metagenomic, Phase:  Gene Finding\n");
                else
                    fprintf(stderr, "Request:  Single Genome, Phase:  Gene Finding\n");
            }

            /* Read and process each sequence in the file in succession */
            sprintf(cur_header, "Prodigal_Seq_1");
            sprintf(new_header, "Prodigal_Seq_2");
            while ((slen = next_seq_multi(input_ptr, seq, useq, &num_seq, &gc,
                                          do_mask, mlist, &nmask, cur_header, new_header)) != -1)
            {
                rcom_seq(seq, rseq, useq, slen);
                if (slen == 0)
                {
                    fprintf(stderr, "\nSequence read failed (file must be Fasta, ");
                    fprintf(stderr, "Genbank, or EMBL format).\n\n");
                    exit(14);
                }

                if (quiet == 0)
                {
                    fprintf(stderr, "Finding genes in sequence #%d (%d bp)...", num_seq, slen);
                }

                /* Reallocate memory if this is the biggest sequence we've seen */
                if (slen > max_slen && slen > STT_NOD * 8)
                {
                    nodes = (struct _node *)realloc(nodes, (int)(slen / 8) * sizeof(struct _node));
                    if (nodes == NULL)
                    {
                        fprintf(stderr, "Realloc failed on nodes\n\n");
                        exit(11);
                    }
                    max_slen = slen;
                }

                /* Calculate short header for this sequence */
                calc_short_header(cur_header, short_header, num_seq);

                if (is_meta == 0)
                { /* Single Genome Version */

                    /***********************************************************************
                      Find all the potential starts and stops, sort them, and create a
                      comprehensive list of nodes for dynamic programming.
                    ***********************************************************************/
                    nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
                    qsort(nodes, nn, sizeof(struct _node), &compare_nodes);

                    /***********************************************************************
                      Second dynamic programming, using the dicodon statistics as the
                      scoring function.
                    ***********************************************************************/
                    score_nodes(seq, rseq, slen, nodes, nn, &tinf, closed, is_meta);
                    if (start_ptr != stdout)
                        write_start_file(start_ptr, nodes, nn, &tinf, num_seq, slen, 0, NULL,
                                         VERSION, cur_header);
                    record_overlapping_starts(nodes, nn, &tinf, 1);
                    ipath = dprog(nodes, nn, &tinf, 1);
                    eliminate_bad_genes(nodes, ipath, &tinf);
                    ng = add_genes(genes, nodes, ipath);
                    tweak_final_starts(genes, ng, nodes, nn, &tinf);
                    record_gene_data(genes, ng, nodes, &tinf, num_seq);
                    if (quiet == 0)
                    {
                        fprintf(stderr, "done!\n");
                    }

                    /* Output the genes */
                    print_genes(output_ptr, genes, ng, nodes, slen, output, num_seq, 0, NULL,
                                &tinf, cur_header, short_header, VERSION);
                    fflush(output_ptr);
                    if (trans_ptr != stdout)
                        write_translations(trans_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                                           &tinf, num_seq, short_header);
                    if (nuc_ptr != stdout)
                        write_nucleotide_seqs(nuc_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                                              &tinf, num_seq, short_header);

                }

                else
                { /* Metagenomic Version */

                    low = 0.88495 * gc - 0.0102337;
                    if (low > 0.65)
                        low = 0.65;
                    high = 0.86596 * gc + .1131991;
                    if (high < 0.35)
                        high = 0.35;

                    max_score = -100.0;
                    for (int i = 0; i < NUM_META; i++)
                    {
                        if (i == 0 || meta[i].tinf->trans_table !=
                                          meta[i - 1].tinf->trans_table)
                        {
                            memset(nodes, 0, nn * sizeof(struct _node));
                            nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask,
                                           meta[i].tinf);
                            qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
                        }
                        if (meta[i].tinf->gc < low || meta[i].tinf->gc > high)
                            continue;
                        reset_node_scores(nodes, nn);
                        score_nodes(seq, rseq, slen, nodes, nn, meta[i].tinf, closed, is_meta);
                        record_overlapping_starts(nodes, nn, meta[i].tinf, 1);
                        ipath = dprog(nodes, nn, meta[i].tinf, 1);

                        if (ipath >= 0 && ipath < nn)
                        {
                            if (nodes[ipath].score > max_score)
                            {
                                max_phase = i;
                                max_score = nodes[ipath].score;
                                eliminate_bad_genes(nodes, ipath, meta[i].tinf);
                                ng = add_genes(genes, nodes, ipath);
                                tweak_final_starts(genes, ng, nodes, nn, meta[i].tinf);
                                record_gene_data(genes, ng, nodes, meta[i].tinf, num_seq);
                            }
                        }
                        // else
                        // {
                        //     fprintf(stderr, "Error: ipath (%d) is out of bounds (nn=%d).\n", ipath, nn);
                        //     exit(15);
                        // }
                    }

                    /* Recover the nodes for the best of the runs */
                    memset(nodes, 0, nn * sizeof(struct _node));
                    nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask,
                                   meta[max_phase].tinf);
                    qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
                    score_nodes(seq, rseq, slen, nodes, nn, meta[max_phase].tinf, closed,
                                is_meta);
                    if (start_ptr != stdout)
                        write_start_file(start_ptr, nodes, nn, meta[max_phase].tinf,
                                         num_seq, slen, 1, meta[max_phase].desc, VERSION,
                                         cur_header);

                    if (quiet == 0)
                    {
                        fprintf(stderr, "done!\n");
                    }

                    /* Output the genes */
                    print_genes(output_ptr, genes, ng, nodes, slen, output, num_seq, 1,
                                meta[max_phase].desc, meta[max_phase].tinf, cur_header,
                                short_header, VERSION);
                    fflush(output_ptr);
                    if (trans_ptr != stdout)
                        write_translations(trans_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                                           meta[max_phase].tinf, num_seq, short_header);
                    if (nuc_ptr != stdout)
                        write_nucleotide_seqs(nuc_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                                              meta[max_phase].tinf, num_seq, short_header);
                }

                /* Reset all the sequence/dynamic programming variables */
                memset(seq, 0, (slen / 4 + 1) * sizeof(unsigned char));
                memset(rseq, 0, (slen / 4 + 1) * sizeof(unsigned char));
                memset(useq, 0, (slen / 8 + 1) * sizeof(unsigned char));
                memset(nodes, 0, nn * sizeof(struct _node));
                nn = 0;
                slen = 0;
                ipath = 0;
                nmask = 0;
                strcpy(cur_header, new_header);
                sprintf(new_header, "Prodigal_Seq_%d\n", num_seq + 1);
            }

            if (num_seq == 0)
            {
                fprintf(stderr, "\nError:  no input sequences to analyze.\n\n");
                exit(18);
            }

            /* Free all memory */
            free(seq);
            free(rseq);
            free(useq);
            free(nodes);
            free(genes);
            for (i = 0; i < NUM_META; i++)
                free(meta[i].tinf);

            /* Close all the filehandles and exit */
            INPUT_CLOSE(input_ptr);
            if (output_ptr != stdout)
                fclose(output_ptr);
            if (start_ptr != stdout)
                fclose(start_ptr);
            if (trans_ptr != stdout)
                fclose(trans_ptr);

            /* Remove tmp file */
            if (piped == 1 && remove(input_copy) != 0)
            {
                fprintf(stderr, "Could not delete tmp file %s.\n", input_copy);
                exit(18);
            }
            if (is_meta == 0){
                cout << "Single Genome: " << fnaFiles[a] << " processed successfully.\n";
            } else {
                cout << "Metagenomic: " << fnaFiles[a] << " processed successfully.\n";
            }
            fileProcessed[a] = true;
            if (std::all_of(fileProcessed.begin(), fileProcessed.end(), [](bool processed)
                            { return processed; }))
            {
                cout << "finished!!";
            }

        }
        }
    if (std::all_of(fileProcessed.begin(), fileProcessed.end(), [](bool processed)
                    { return processed; }))
    {
        return 0;
    }

    return 0;
}
size_t getSequenceLength(fptr fp)
{
    char line[MAX_LINE + 1];
    size_t totalLength = 0;

    line[MAX_LINE] = '\0';
    while (INPUT_GETS(line, MAX_LINE, fp) != NULL)
    {
        if (line[0] == '>')
            continue;

        for (unsigned int i = 0; i < strlen(line); i++)
        {
            if ((line[i] >= 'A' && line[i] <= 'Z') || (line[i] >= 'a' && line[i] <= 'z'))
            {
                // cout << line[i] << flush;
                totalLength++;
            }
        }
    }
    cout << "Total sequence length: " << totalLength << endl;
    return totalLength;
}

bool hasReadPermission(const std::string &filePath)
{
    struct stat fileStat;
    if (stat(filePath.c_str(), &fileStat) != 0)
    {
        std::cerr << "Error: Could not access file: " << filePath << "\n";
        return false;
    }

    // Check if the file has read permission
    if ((fileStat.st_mode & S_IRUSR) || (fileStat.st_mode & S_IRGRP) || (fileStat.st_mode & S_IROTH))
    {
        return true;
    }
    else
    {
        std::cerr << "Error: File does not have read permission: " << filePath << "\n";
        return false;
    }
}
