#!/usr/bin/env python

import sys
import os
import argparse

# What does this script do? 

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.abspath(os.path.expanduser(values))]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_arguments(): #edit as needed!
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Appends vcf with antifungal gene info")
    parser.add_argument('-v', '--vcf',
        help ='VCF file or any file with chrom and pos as\
        the first two columns', 
        type=is_file)
    return parser.parse_args()

args = get_arguments()

autosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

def make_antifungaldict(autosomes):
    antifungals = {}
    for chromosome in autosomes:
        antifungals[chromosome] = {}
    with open("AFimmunity_genes_unified_coordinates.txt", 'r') as targets:
        for line in targets:
            if line[0] != "#":
                line = line.strip().split()
                chrom,start,stop,name = line[0],int(line[1]),int(line[2]),line[3]
                if chrom in autosomes:
                    antifungals[chrom][name] = [start,stop,name]
    return(antifungals)

def make_AFregDict(autosomes):
    AFregDict = {}
    for chromosome in autosomes:
        AFregDict[chromosome] = {}
    with open("AFimmunity_reg_regions.txt", 'r') as REG:
        for line in REG:
            if line[0] != "#":
                line = line.strip().split()
                chrom,start,stop,name = line[0],int(line[1]),int(line[2]),line[3]
                if chrom in autosomes:
                    AFregDict[chrom][name] = [start,stop,name]
    return(AFregDict)

def append(autosomes, antifungals, AFregDict):
    outfilename = str(args.vcf).split(".")[0] + '_antifungal.txt'
    with open(args.vcf, "r") as vcf, open(outfilename, 'w') as outfile:
        for line in vcf:
            if line[0] != "#":
                line = line.strip().split('\t')
                chrom = "chr" + line[0]
                bp = line[1]
                if chrom in autosomes:
                    candReg = []
                    candGene = []
                    for v in antifungals[chrom].values():
                        if int(bp) >= int(v[0]) and int(bp) <= int(v[1]):
                            candGene.append(v[2])
                    for val in AFregDict[chrom].values():
                        if int(bp) >= int(val[0]) and int(bp) <= int(val[1]):
                            candReg.append(val[2])
                    if not candReg:
                        candReg = ["NA"]
                    if not candGene:
                        candGene = ["NA"]
                    outfile.write("\t".join(line) + "\t" + ",".join(candReg) + '\t' + ",".join(candGene) + "\n")

antifungals = make_antifungaldict(autosomes)
AFregDict = make_AFregDict(autosomes)
append(autosomes, antifungals, AFregDict)
