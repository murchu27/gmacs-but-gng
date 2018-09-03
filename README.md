# gmacs-but-gng
A new implementation of GMACS which aims to cluster using GNG maps, rather than SOM maps. 

_Original work on GMACS is credited to Pilib O'Broin (pilib.obroin@nuigalway.ie)_

## Example usage

```
g++ gmacs.cpp make_kfv.cpp -o gmacs -std=c++11

./gmacs -f mmset_activated_prom_seqs.fasta
```

GMACS aims to cluster motif files, and sequence files, and naturally a certain format is expected. See outlines below:

### Motif files

- Motif names should be prepended with a `>`.
- Below each name should be the position weight matrix (PWM) for that motif, in which:
  - the rows correspond to the positions of the motif.
  - the columns correspond to the frequencies of the nucleotides A, C, G, T.
- The maximum allowed length of the PWM is set in [gmacs.h], by the value of `MAX_L`.

#### Example (from [jasp79_mat.txt])
```
>AGL3_MADS
0.000	0.969	0.010	0.021
0.031	0.773	0.000	0.196
...(4xl frequency values)
>ARNT_bHLH
0.200	0.800	0.000	0.000
0.950	0.000	0.050	0.000
...
```

### Sequence files

- Sequence names should be prepended with a `>`.
- Below each name should be the actual sequence.
- The length of the sequence is set in [gmacs.h], by the value of `SEQ_L`. 
  - Currently, this MUST be changed if any sequences are of a different length (using the `-s` option when running GMACS). 
- Note that GMACS ignores newline characters when reading sequence files; the only delimiter between sequences is the `>` character.

#### Example (from [jasp79_mat.txt])
```
>chr1_33274558_33277058
AACCTGGGCCATAAAGCTAATTGCCTCCATTCAGCCTGTAATCTAGATGC
CTCGGCCCAGGGCTCACTCTACCCCTCAAACCTTTCCAGTCTCTAGACCC
...(SEQ_L chars per chromosome)
>chr1_150232838_150235338
GGGCGGAGCTCACCTTGGCCGAGGCGCGGCGGACGCTGGGCGAGCTGGGC
...
```
