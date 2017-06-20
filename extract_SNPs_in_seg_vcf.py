#!/usr/bin/env python

import argparse
import sys
import os
from collections import defaultdict

#########################################################################
###This script was written some time ago and is designed to extract SNPs
###identified as being in consensus ROH by the PLINK algorithm and write
###them to file. Updated 6.20.2017 with argparse.
#########################################################################

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Extract SNPs in consensus ROH")
    parser.add_argument('-c', '--consensusROH', required=True,
        help = 'Consensus ROH file (generally .overlap)',
        type = is_file)
    parser.add_argument('-v', '--vcf', required=True,
        help = 'VCF',
        type = is_file)
    parser.add_argument('-o', '--output', required=True,
        help = 'outfile name')
    return parser.parse_args()

args = get_arguments()

def make_chrdict():
    #make a dictionary of dictionaries, one dictionary per autosome with the geneaccession number as key
    dictionaries = {}
    autosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    for chromosome in autosomes:
        dictionaries[chromosome] = {}
    with open(args.consensusROH, 'r') as CON:
        for line in CON:
            line = line.strip().split()
            seg,startSEG,stopSEG = line[0],str(line[7]),str(line[8])
            chrom = "chr" + line[4]
            if chrom in autosomes:
                dictionaries[chrom][seg] = [startSEG,stopSEG]
    return(dictionaries,autosomes)

def make_gene_segment_overlap(dictionaries,autosomes):
    with open(args.vcf, "r") as VCF, open(args.output, "w") as outfile:
        for line in VCF:
            if line[0] != "#":
                line = line.strip().split()
                chrom,bp = line[0],line[1]
                if chrom in autosomes:
                    for v in dictionaries[chrom].values():
                        if int(bp) >= int(v[0]) and int(bp) <= int(v[1]):
                            outfile.write("\t".join(line) + "\n")           

dictionaries,autosomes = make_chrdict()
for chrom in autosomes:
    print(chrom + "\t" + str(len(dictionaries[chrom])))
make_gene_segment_overlap(dictionaries,autosomes)
