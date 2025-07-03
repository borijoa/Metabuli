#include <sys/stat.h>
#include <unistd.h>
#include "bitmap.h"
#include "dprog.h"
#include "gene.h"
#include "fptr.h"
#include "metagenomic.h"
#include "node.h"
#include "prodigalsequence.h"
#include "training.h"
#include "printBinary.h"
#include <vector>
#include <cstdint>

#define MIN_SINGLE_GENOME 20000
#define IDEAL_SINGLE_GENOME 100000
#define VERSION "2.6.3"
#define DATE "February, 2016"

class Prodigal
{
private:
    int rv, slen, ipath, *gc_frame, output;
    int closed, do_mask, nmask, force_nonsd, user_tt, num_seq, quiet;
    int max_slen, fnum;
    double gc, low, high;
    unsigned char *seq, *rseq, *useq;
    char *train_file, *start_file, *trans_file, *nuc_file;
    char *input_file, *output_file, input_copy[MAX_LINE];
    char cur_header[MAX_LINE], new_header[MAX_LINE], short_header[MAX_LINE];
    FILE *start_ptr, *trans_ptr, *nuc_ptr;
    struct stat fbuf;
    pid_t pid;
    struct _training tinf;
    mask mlist[MAX_MASKS];
    struct _metagenomic_bin *meta;

    /* Set a bit to 1 */
    void set(unsigned char *bm, int ndx)
    {
        bm[ndx >> 3] |= (1 << (ndx & 0x07));
    }

public:
    int fng, ng;
    int is_meta;
    int is_first_meta;
    int nn, max_phase;
    double max_score;
    struct _node *nodes;
    struct _gene *genes;
    struct _gene *finalGenes;
    void ensureDirectoryExists(const std::string &path);
    size_t getSequenceLength(fptr fp);
};
