# primerfind
A python code for demultiplexing. Finds barcodes and tags with error of 2 nt.

Input:
At the moment, reads in barcodes from a default file called "Sample_primers_barcodes.csv" that contains the following columns: Sample, Forward primer, F primer barcode, R primer barcode and Reverse primer. It is essential that any input file has the same structure and naming of the file.

Sequence file is read in from standard input and should be in FASTA format.

Output:

Redirects demultiplexed sample files into a new directory. Outputs also a summary (part of all samples) "results"-file with sample names and numbers of forward and reverse sequences. Furthermore, outputs  a full results file with all the samples and numbers of sequences. 