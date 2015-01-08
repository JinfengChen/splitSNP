splitSNP
========

*splitSNP* is a small utility to extract reads from a BAM file that cover a specified SNP and output the reference and alternate allele containing reads to separate BAM files.

I created *splitSNP* for the analysis of allele-specific methylation from bisulfite amplicon experiments sequenced on the MiSeq.

Usage
---------
    splitSNP.py [-h] [--max_depth MAX_DEPTH] input_bam output_prefix SNP


* input_bam - Input BAM file
* output_prefix - Prefix for output files
* SNP - SNP position, reference and alternate allele of interest in chr:position:ref:alt format, eg chr21:11106932:A:G. For deletion analysis the ref should be 'D' and alt be the size of the deletion in basepairs, eg chr11:67351213:D:64 .Coordinates are 1-based.
* --max_depth MAX_DEPTH -  Maximum number of reads to process at the specified SNP position (defaults to 1000000)


SNP usage
---------
This example uses *splitSNP* to extract and separate the **A** (reference allele) and **G** (alternate allele) containing reads spanning the SNP *rs11184058* which located at chr21, position 11,106,932 in hg19 (obtained from UCSC genome browser).

    splitSNP.py MiSeq_20.bam MiSeq_20.rs11184058 chr21:11106932:A:G

<!-- separate input and output -->

    BAM file 'MiSeq_20.bam' does not have an index, creating one...
    Processing 118528 reads covering SNP position chr21:11106932 in MiSeq_20.bam
    108970 reads with the reference allele 'A' written to MiSeq_20.rs11184058.ref.A.bam
    9254 reads with the alternate allele 'G' written to MiSeq_20.rs11184058.alt.G.bam
    Discarded 237 'C' alleles
    Discarded 42 'T' alleles
    Discarded 35 'N' alleles

Deletion usage
--------------
This example uses *splitSNP* to extract reads from the human *GSTP1* locus, separating reads matching the reference sequence from those containing an engineered 64 basepair deletion spanning chr11:67,351,213-67,351,276 (hg19).



    splitSNP.py MiSeq_1_S1.bam MiSeq_1_S1.GSTP1 chr11:67351213:D:64

<!-- separate input and output -->

    7120 reads with the reference sequence written to MiSeq_1_S1.GSTP1.ref.bam
    2466 reads with the 64bp deletion written to MiSeq_1_S1.GSTP1.del.64bp.bam


TODOs
-----
* There's nothing here!
