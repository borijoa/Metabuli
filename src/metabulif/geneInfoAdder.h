#ifndef DERIVED_INDEXCREATOR_H
#define DERIVED_INDEXCREATOR_H

#include "IndexCreator.h"

class geneInfoAdder : public IndexCreator {
protected:
    //Inputs
    std::string unirefId2IdxFileName;
    std::string unirefIdx2taxIdFileName;
    std::string ncbi2gtdbFileName;

    //Outputs
    std::string geneListFileName;
    size_t fillGeneKmerBuffer(TargetKmerBuffer &kmerBuffer, bool *checker, size_t &processedSplitCnt,
    const LocalParameters &par);
    
    //Required to annotate with Uniref90
    std::vector<string> uniRefIds;
    std::vector<TaxID> uniRefTaxIds;
    std::unordered_set<TaxID> speciesTaxIds;
    std::unordered_map<TaxID, TaxID> ncbi2gtdb;

    // Results
    std::unordered_map<uint32_t, TaxID> regionId2taxId;
    std::unordered_map<uint32_t, uint32_t> regionId2unirefId;

public:
    geneInfoAdder(const LocalParameters & par, TaxonomyWrapper * taxonomy)
        : IndexCreator(par, taxonomy) {}

    ~geneInfoAdder();

    // search from Uniref90
    void customMethod() {
        std::cout << "This is a custom method in DerivedIndexCreator." << std::endl;
    }

    // 
    void createGeneIndex(const LocalParameters & par);
    
};

#endif // DERIVED_INDEXCREATOR_H
