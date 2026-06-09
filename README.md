# The Extended Nei-Gojobori Method for Calculating Synonymous and Nonsynonymous Distances

## Introduction

This program is designed for the incorporating mutation spectrum into the Nei-Gojobori method (Nei and Gojobori, 1986) by using the relative rates of six mutation types (AT -> TA, AT -> CG, AT -> GC, CG -> TA, CG -> AT, CG -> GC) to calculate S and N for a sequence. 

This program is written in Python and requires Python >= 3.10.

The program should be cited with:

<citation to be added>


## Usage

The program can be run from the command line:

`python </path/run_extended_ng.py> --sequence-file <path> --outfile <path> --optional-argument <optional_argument-value>`

Required arguments:

`--sequence-file <path>`

`--outfile <path>`

Optional arguments:

`--spectrum-file <path>`

Specify mutation spectrum to be used.
	
If not specified, the program assumes all six mutation types to have equal mutation rates, as did in Nei and Gojobori, 1986.

`--genetic-code <int>`

Specify genetic code to be used. 
	
Please use the index of genetic code as listed in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.
	If not specified, the standard code (--genetic-code 1) is assumed.

`--errorfile <path>`

If specified, errors and warnings will be printed into the error file.
	
If not, errors and warnings will appear in terminal.


## Sequence file

The sequence file should consist a header line containing two integers:

<number_of_sequences> <sequence_length>

The remaining lines should contain one sequence name followed by one sequence line. Therefore, each sequence is represented by two lines:

<sequence_name>
	
<sequence>

Requirements for the sequence file:

1. The sequence length must be a multiple of 3.
2. All sequences must have the same length as specified in the first line.
3. The number of sequences must match the number specified in the first line.
4. Sequences should be coding DNA sequences and should be aligned by codon.
5. Sequences should only contain A, T, C, G and a, t, c, g.
6. Stop codons should be removed from sequences.

See rnase.seq for an example.


## Spectrum file

The spectrum file has four possible formats.

### Format 1: ordered 6-parameter format

There should be 6 lines in the file. Each line should be a float describing:

	ATtoTA relative mutation rate
    ATtoCG relative mutation rate
    ATtoGC relative mutation rate
    CGtoTA relative mutation rate
    CGtoAT relative mutation rate
    CGtoGC relative mutation rate

An example: 

    0.33
    0.33
    0.34 
    0.20
    0.50
    0.30

These numbers could be actual mutation rate times some arbitrary constant, as long as these 6 parameters are on the same scale.

### Format 2: labeled 6-parameter format

An example:

    ATtoTA 0.33
    CGtoGC 0.30
    ATtoCG 0.33
    CGtoAT 0.50
    ATtoGC 0.34
    CGtoTA 0.20

Please note that labels should be consistent with indicated in the example. The order of labels does not matter.

These numbers could be actual mutation rate times some arbitrary constant, as long as these 6 parameters are on the same scale.

### Format 3: ordered 8-parameter format

There should be 8 lines in the file. Each line should be a float describing:

	ATtoTA relative mutation rate (a)
    ATtoCG relative mutation rate (b)
    ATtoGC relative mutation rate (c)
    CGtoTA relative mutation rate (d)
    CGtoAT relative mutation rate (e)
    CGtoGC relative mutation rate (f)
	total A/T mutation rate (g)
	total C/G mutation rate (h)

An example:

    1/3
    1/3
    1/3
    1/2
    1/4
    1/4
    3
    4.5

These numbers could be actual mutation rate times some arbitrary constant, as long as:

 a, b, and c should be on the same scale. 
 
 d, e, and f should be on the same scale.
 
 g and h should be on the same scale.

### Format 4: labeled 8-parameter format

An example:

    ATtoTA 1
    ATtoCG 1
    ATtoGC 1
    CGtoTA 0.4
    CGtoAT 0.5
    CGtoGC 0.3
    ATrate 0.8
    CGrate 1.2

Please note that labels should be consistent with indicated in the example. The order of labels does not matter.

These numbers could be actual mutation rate times some arbitrary constant, as long as:

 Relative mutation rate of different point mutations from A/T should be on the same scale. 
 
 Relative mutation rate of different point mutations from C/G should be on the same scale.
 
 Mutation rate of A/T and C/G should be on the same scale.


## Outfile 

The outfile may include the following section:
1. Genetic code
2. Number of taxa 
3. Sequence names 
4. Information about removed codons 
5. Number of nonsynonymous differences between two sequences (n) 
6. Number of synonymous differences between two sequences (s) 
7. Mutation spectrum
8. Number of nonsynonymous sites (N) and synonymous sites (S) of a sequence
9. Nonsynonymous p-distances between two sequences (pN) 
10. Synonymous p-distances between two sequences (pS)
11. pN/pS
12. Jukes-Cantor corrected nonsynonymous distances between two sequences (dN) 
13. Jukes-Cantor corrected synonymous distances between two sequences (dS)
14. dN/dS

See rnase.out for an example.
