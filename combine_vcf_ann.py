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
        help = 'SNPs.overlap - vcf file with snps of interest',
        type=is_file)
    parser.add_argument('-tg', '--tgFile',
        help ='Thousand genomes vcf', 
        type=is_file)
    return parser.parse_args()

args = get_arguments()

def make_snpDict():
    """Make a dictionary with the chr_pos as the key and 
    an empty list as value"""
    snpDict = {}
    with open(args.tgFile, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            chrom = line[0]
            pos = line[1]
            snpID = chrom + '_' + pos
            snpDict[snpID] = [line[2], line[7]]
    return snpDict

def make_snpDictAlt():
    """Make a dictionary with the chr_pos as the key and 
    an empty list as value"""
    snpDict = {}
    cols = ['AF','EUR_AF','EAS_AF','SAS_AF','AFR_AF','AMR_AF','AA','tri','CSQ','FUNSEQ']
    with open(args.tgFile, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            chrom = line[0]
            pos = line[1]
            ALTs = line[4].split(",")
            INFOs = line[7].split(";")
            for i,allele in enumerate(ALTs):
                snpID = chrom + '_' + pos + '_' + ALTs[i]
                snpDict[snpID] = {}
                snpDict[snpID]['rs'] = line[2]
                snpDict[snpID]['info'] = line[7]
                for x in cols:
                    snpDict[snpID][x] = "NA"
                for field in INFOs:
                    try:
                        header,value = field.split("=")
                        values = value.split(",")
                        if len(values) > 1:
                            snpDict[snpID]['tri'] = "TRUE"
                            if header == 'CSQ':
                                for x,y in enumerate(values):
                                    if y[0] == ALTs[i]:
                                        val = values[x]
                            else:
                                val = values[i]
                        else:
                            snpDict[snpID]['tri'] = "FALSE"
                            val = value
                        snpDict[snpID][header] = val
                    except ValueError:
                        print field
    return snpDict

def write_file(snpDict):
    outfileName = args.snpFile.split(".")[0] + '_mod.vcf'
    with open(args.snpFile, 'r') as infile, open(outfileName, 'w') as outfile:
        for line in infile:
            line = line.strip().split('\t')
            chrom = line[0].strip("chr")
            pos = line[1]
            alt = line[4]
            triTest = alt.split(",")
            if len(triTest) > 1:
                print("{0} {1} has more than one alt allele".format(chrom, pos))
            snpID = chrom + '_' + pos 
            if snpID in snpDict:
                if snpDict[snpID][0] == line[2]:
                    outfile.write("\t".join(line[0:7]) + '\t' + snpDict[snpID][1] + '\t' + "\t".join(line[8:]) + '\n')
                elif snpDict[snpID][0][0:2] == "rs" and line[2] == ".":
                    outfile.write("\t".join(line[0:2]) + '\t' + snpDict[snpID][0] + '\t' + "\t".join(line[3:7]) + '\t' + snpDict[snpID][1] + '\t' + "\t".join(line[8:]) + '\n')
                else:
                    print('Disagreement between rs#s at {0}'.format(snpID))
                    outfile.write("\t".join(line[0:2]) + '\t' + line[2] + ',' + snpDict[snpID][0] + '\t' + "\t".join(line[3:7]) + '\t' + snpDict[snpID][1] + '\t' + "\t".join(line[8:]) + '\n')
            else:
                outfile.write("\t".join(line) + '\n')

def write_file_alt(snpDict):
    outfileName = args.snpFile.split(".")[0] + '_mod.vcf'
    with open(args.snpFile, 'r') as infile, open(outfileName, 'w') as outfile:
        for line in infile:
            line = line.strip().split('\t')
            chrom = line[0].strip("chr")
            pos = line[1]
            alt = line[4]
            triTest = alt.split(",")
            if len(triTest) > 1:
                print("{0} {1} has more than one alt allele".format(chrom, pos))
            snpID = chrom + '_' + pos + '_' + alt
            if snpID in snpDict:
                if snpDict[snpID]['rs'] == line[2]:
                    outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                        ("\t".join(line[0:7]),
                        snpDict[snpID]['info'], 
                        "\t".join(line[8:]), 
                        snpDict[snpID]['AF'],
                        snpDict[snpID]['EUR_AF'],
                        snpDict[snpID]['EAS_AF'],
                        snpDict[snpID]['SAS_AF'],
                        snpDict[snpID]['AFR_AF'],
                        snpDict[snpID]['AMR_AF'],
                        snpDict[snpID]['AA'],
                        snpDict[snpID]['tri'],
                        snpDict[snpID]['CSQ'],
                        snpDict[snpID]['FUNSEQ']))
                elif snpDict[snpID]['rs'][0:2] == "rs" and line[2] == ".":
                    outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                        ("\t".join(line[0:2]),
                        snpDict[snpID]['rs'],
                        "\t".join(line[3:7]),
                        snpDict[snpID]['info'], 
                        "\t".join(line[8:]), 
                        snpDict[snpID]['AF'],
                        snpDict[snpID]['EUR_AF'],
                        snpDict[snpID]['EAS_AF'],
                        snpDict[snpID]['SAS_AF'],
                        snpDict[snpID]['AFR_AF'],
                        snpDict[snpID]['AMR_AF'],
                        snpDict[snpID]['AA'],
                        snpDict[snpID]['tri'],
                        snpDict[snpID]['CSQ'],
                        snpDict[snpID]['FUNSEQ']))
                else:
                    print('Disagreement between rs#s at {0}'.format(snpID))
                    outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                        ("\t".join(line[0:2]),
                        line[2] + ',' + snpDict[snpID]['rs'],
                        "\t".join(line[3:7]),
                        snpDict[snpID]['info'], 
                        "\t".join(line[8:]), 
                        snpDict[snpID]['AF'],
                        snpDict[snpID]['EUR_AF'],
                        snpDict[snpID]['EAS_AF'],
                        snpDict[snpID]['SAS_AF'],
                        snpDict[snpID]['AFR_AF'],
                        snpDict[snpID]['AMR_AF'],
                        snpDict[snpID]['AA'],
                        snpDict[snpID]['tri'],
                        snpDict[snpID]['CSQ'],
                        snpDict[snpID]['FUNSEQ']))
            else:
                outfile.write("\t".join(line) + 'NA\t'*10 + 'NA\n')

#snpDict = make_snpDict()
#write_file(snpDict)
snpDict = make_snpDictAlt()
write_file_alt(snpDict)
