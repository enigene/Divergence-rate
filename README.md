# Divergence rate
Calculate pairwise sequence divergence rate beween aligned FASTA sequnces

The program requires an alignment in FASTA format.

### Usage:
`awk -v rules=ggsnns -f divergr.awk alignment.fas > output.txt`

### Rules:
+ gg — gap-gap reduce length;
+ gs — gap-site reduce length;
+ nn — N-N not counted as similarity;
+ ns — N-site reduce length.

Multiple rules can be used:
+ ggs — gap-gap and gap-site reduce length;
+ ggsnns includes all of the above.

### Options:
+ `-v printDivAll=1` — prints divergence for each sequence with all sequences, except itself (default);
+ `-v printSeq=1` — prints mean and median divergence for each sequences;
+ `-v printTotal=1` — prints minimal, maximum, mean, and median divergence for all sequences;
+ `-v debug=1` — prints variables, arrays and statements during execution.
