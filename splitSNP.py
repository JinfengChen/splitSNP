#! /usr/bin/env python
def main():
    import argparse
    import pysam
    import os
    import tempfile
    from collections import Counter

    def can_create_file(folder_path):
        try:
            tempfile.TemporaryFile(dir=folder_path)
            return True
        except OSError:
            return False

    parser = argparse.ArgumentParser(description="splitSNP.py - "
        "Extracts reads from a BAM file that cover a specified SNP and "
        "writes the reference and alternate allele containing reads to separate BAM files.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("SNP", help="SNP position, reference and alternate "
        "allele of interest in chr:position:ref:alt format, eg chr21:11106932:A:G. "
        "Coordinates are 1-based.")
    parser.add_argument("--max_depth", type=int, default=1000000, help="Maximum number "
        "of reads to process at the specified SNP position")
    args = parser.parse_args()

    if args.max_depth < 8000:
        print("Specified max_depth is too low - changing to 8000")
        args.max_depth = 8000

    # Check input file exists, and thet output folder is writeable
    if not os.path.isfile(args.input_bam):
        print("Input BAM file %s does not exist!" % args.input_bam)
        return

    if not can_create_file(os.path.dirname(args.output_prefix)):
        print("Output path %s is not writable!" % os.path.dirname(args.output_prefix))
        return

    # Check alleles are valid
    nucleotides = ['A', 'C', 'G', 'T']
    try:
        chrom, pos, ref, alt = args.SNP.split(":")
    except ValueError:
        print("SNP specified '%s' not in chr:position:ref:alt format" % args.SNP)
        return

    ref = ref.upper()
    alt = alt.upper()
    if ref not in nucleotides:
        print("Reference allele %s is not A, C, G or T" % ref)
        return
    if alt not in nucleotides:
        print("Alternate allele %s is not A, C, G or T" % alt)
        return

    # Check validity of chrom
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    chrom_no = next((i for i in range(len(samfile.references)) if samfile.references[i]==chrom), None)
    if chrom_no is None:
        print("Chromosome '%s' not in BAM '%s'" % (chrom, args.input_bam))
        return
    
    # Check validity of pos
    try:
        pos = int(pos)
    except ValueError:
        print("Position '%s' is not valid" % pos)
        return
    if pos >= samfile.lengths[chrom_no]:
        print("Position '%s' is out of bounds of chromosome '%s'" % (pos, chrom))
        return

    # Index samfile if one does not already exist
    if not samfile._hasIndex():
        print("BAM file '%s' does not have an index, creating one..." % args.input_bam)
        samfile.close()
        pysam.index(args.input_bam)
        samfile = pysam.AlignmentFile(args.input_bam, "rb")

    # PASS 1 - find readnames of all reads with ref/alt alleles
    alleles = list()
    ref_readnames = set()
    alt_readnames = set()
    for pileup in samfile.pileup(chrom, pos-1, pos, max_depth=args.max_depth):
        if pileup.reference_pos == pos-1: # filter for position of interest
            print "Processing %s reads covering SNP position %s:%s in %s" % (
                len(pileup.pileups), chrom, pos, args.input_bam)
            for read in pileup.pileups:
                SNP_base = read.alignment.query_sequence[read.query_position]
                alleles.append(SNP_base)
                if SNP_base == ref:
                    ref_readnames.add(read.alignment.query_name)
                elif SNP_base == alt:
                    alt_readnames.add(read.alignment.query_name)

    # Remove reads in both
    ref_and_alt = ref_readnames.intersection(alt_readnames)
    if len(ref_and_alt) > 0:
        ref_readnames = ref_readnames.difference(ref_and_alt)
        alt_readnames = alt_readnames.difference(ref_and_alt)

    # PASS 2 - output reads matching above readnames to two new bamfiles
    ref_filename = args.output_prefix+".ref."+ref+".bam"
    alt_filename = args.output_prefix+".alt."+alt+".bam"
    ref_bam = pysam.AlignmentFile(ref_filename, "wb", template=samfile)
    alt_bam = pysam.AlignmentFile(alt_filename, "wb", template=samfile)
    ref_count = 0
    alt_count = 0
    for read in samfile.fetch(chrom, pos-1, pos):
        if read.query_name in ref_readnames:
            ref_bam.write(read)
            ref_count += 1
        elif read.query_name in alt_readnames:
            alt_bam.write(read)
            alt_count += 1

    ref_bam.close()
    alt_bam.close()
    samfile.close()

    # Print a summary
    allele_count = Counter(alleles)
    print "%s reads with the reference allele '%s' written to %s" % (
        ref_count, ref, ref_filename)
    print "%s reads with the alternate allele '%s' written to %s" % (
        alt_count, alt, alt_filename)
    for x in allele_count:
        if x != ref and x != alt:
            print "Discarded %s '%s' alleles" % (allele_count[x], x)

if __name__ == '__main__':
    main()
