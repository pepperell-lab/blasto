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
        description="Annotes CADD or modCADD file with GTEx eQTLs")
    parser.add_argument('-i', '--snpFile',
        help ='CADD annotation file or any file with chrom and pos as\
        the first two columns', 
        type=is_file)
    parser.add_argument('-d', '--directory',
        help ='directory containing GTEx eQTL files', 
        type=is_dir)
    return parser.parse_args()

args = get_arguments()

def make_snpDict(snpIN):
    """Make a dictionary with the chr_pos as the key and 
    an empty list as value"""
    snpDict = {}
    with open(snpIN, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            chrom = line[0]
            pos = line[1]
            snpID = chrom + '_' + pos
            snpDict[snpID] = []
    return snpDict

def add_GTEx(snpDict, dirIN):
    for root, dirs, files in os.walk(dirIN):
        for tissueFile in files:
            tissue = "".join(tissueFile.split("_")[:-1])
            with open(os.path.join(root,tissueFile), 'r') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    chrom = line[13]
                    pos = line[14]
                    snpID = chrom + '_' + pos
                    geneName = line[26]
                    if snpID in snpDict:
                        eQTL = tissue + '-' + geneName
                        snpDict[snpID].append(eQTL)
    return snpDict
        
def write_file(snpDict, snpIN):
    outfileName = snpIN.split(".")[0] + '_GTEx.txt'
    with open(args.snpFile, 'r') as infile, open(outfileName, 'w') as outfile:
        for line in infile:
            line = line.strip().split('\t')
            chrom = line[0]
            pos = line[1]
            snpID = chrom + '_' + pos
            if len(snpDict[snpID]) > 0:
                outfile.write("\t".join(line) + "\t" + ",".join(snpDict[snpID]) + "\n")
            else:
                outfile.write("\t".join(line) + "\t" + "NA" + "\n")

snpDict = make_snpDict(args.snpFile)
snpDict = add_GTEx(snpDict, args.directory)
write_file(snpDict, args.snpFile)
