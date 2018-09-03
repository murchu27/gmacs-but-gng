//include guard
#ifndef GMACS_H
#define GMACS_H

//include directives
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <cctype>

//definitions of kfv forming/clustering parameter values
#define POP 100
#define GEN 300
#define SP 1.00
#define MUTATE 0.05
#define TRIM 0.30
#define REP 0.05
#define KMER 4
#define LMER 10
#define MIN_MOTIF 4
#define MAX_L 30
#define SEQ_L 2500
#define FN 100

//motif struct to store motifs being read in
struct motif
{
    char name[FN];
    int length;
    double pssm[MAX_L][4];
    double ic[MAX_L];
};

//sequence struct to store sequences being read in
struct sequence
{
    char name[FN];
    char seq[SEQ_L];
};

//declare function/utility methods
void read_file(char* file);
void mot_read(FILE* mot);
void seq_read(FILE* seq);
void trim_motif(struct motif*);
void get_ic(struct motif* motif_x);
void mot_to_kfv(motif* motifs);
void seq_to_kfv(sequence* sequences);
void init_kfv();
void report_kfv();
void free_kfv();

//declare parameters for reading/writing motifs/sequences
extern FILE* out;
extern FILE* outR;
extern int DN;
extern char file[FN];

//declare parameters for kfv forming & clustering
extern double **kfv_f, **kfv_r;
extern int pop, gen;
extern double mut, trim, rep;
extern int trim_i, mut_i;
extern int min_motif;
extern int kmer, lmer, seql;
extern int kfv_size;

#endif //GMACS_H
