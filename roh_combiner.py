#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt
import glob
import subprocess
import shlex

####################################################################
##This script is written to unify ROH identified by bcftools and 
##garlic.
####################################################################


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

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

def get_arguments(): 
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Combine ROH identified by Garlic and Bcftools")
    parser.add_argument('-g', '--garlic', required=True,
        help = 'Garlic ROH output *roh.bed', 
        type = is_file)
    parser.add_argument('-b', '--bcftools', required=True,
        help = 'Bcftools ROH output - per region - RG',
        type = is_file)
    parser.add_argument('-o', '--outBCF', required=True,
        help = 'Combined ROH')
    parser.add_argument('-o2', '--outGAR', required=True,
        help = 'Combined ROH')
    return parser.parse_args()

args = get_arguments()

autosomes = []
for i in range(1,23):
        chrom = "chr" + str(i)
        autosomes.append(chrom)
        
individuals = ["BDH001", "BDH002", "BDH003", "BDH004", "BDH005", "BDH006", "BDH007", "BDH008", "BDH009", "BDH010", "BDH018"]

def add_bcftools_roh():
    with open(args.bcftools, 'r') as infile, open(args.outBCF, 'w') as outfile:
        next(infile)
        for line in infile:
            line = line.strip().split('\t')
            samp = line[1].replace("-","")
            chrom = "chr" + line[2]
            start = int(line[3])
            stop = int(line[4])
            if chrom in autosomes:
                outfile.write( 
                chrom + '\t' +
                str(start) + '\t' +
                str(stop) + '\t' +
                samp + '\t' +
                "bcftools" + '\t' +
                line[5] + '\t' +
                "NA" + '\t' +
                line[6] + '\t' +
                line[7] + '\t' +
                "NA\tNA\tNA\n")

def add_garlic_roh():
    samp = None
    with open(args.garlic, 'r') as infile, open(args.outGAR, 'w') as outfile:
        for line in infile:
            line = line.strip().split('\t')
            if len(line) == 1:
                header = "".join(line).split()
                samp = header[2]
                print(samp)
            elif len(line) == 9:
                chrom = line[0]
                start = int(line[1])
                stop = int(line[2])
                outfile.write(
                chrom + '\t' + 
                str(start) + '\t' +
                str(stop) + '\t' +
                samp + '\t' +
                "garlic" + '\t' +
                line[4] + '\t' +
                line[3] + '\t' +
                "NA\tNA\tNA\tNA\tNA\n")

add_bcftools_roh()
add_garlic_roh()
