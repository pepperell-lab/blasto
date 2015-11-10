#!/usr/bin/env python

from collections import defaultdict

#####################################################################
##This script was written some time ago for exploration of the pilot
##data for the blast project. I have recently (Nov 2015) made a few 
##edits of file paths and outputs to update previous annotations.
##As is, this script serves one purpose: to add annotations to the 
##output from the CADD server. At this point in time I do not plan
##to make it more reusable. Thus, examination of the file paths
##should be done to ensure the desired sources are being drawn from.
####################################################################

autosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

def make_chrdict(autosomes):
    #make a dictionary of dictionaries, one dictionary per autosome with the geneaccession number as key
    dictionaries = {}
    snpIDdict = {}
    for chromosome in autosomes:
        dictionaries[chromosome] = {}
        snpIDdict[chromosome] = {}
    with open("/opt/external_drive/oneill/human/automated/mother_files/pilot_mod.tped", 'r') as genotypes:
        for line in genotypes:
            if line[0] != "#":
                line = line.strip().split()
                chrom = "chr" + line[0]
                loc = line[3]
                snpID = line[1]
                genotypes = line[4:16]
                if chrom in autosomes:
                    dictionaries[chrom][loc] = genotypes
                    snpIDdict[chrom][loc] = snpID
    return(dictionaries,snpIDdict)

def make_antifungaldict(autosomes):
    #make a dictionary of dictionaries, one dictionary per autosome with the geneaccession number as key
    antifungals = {}
    for chromosome in autosomes:
        antifungals[chromosome] = {}
    with open("/opt/external_drive/oneill/antifungal/AFimmunity_genes_unified_coordinates.txt", 'r') as targets:
        for line in targets:
            if line[0] != "#":
                line = line.strip().split()
                chrom,start,stop,name = line[0],int(line[1]),int(line[2]),line[3]
                if chrom in autosomes:
                    antifungals[chrom][name] = [start,stop,name]
    return(antifungals)

def make_AFregDict(autosomes):
    #make a dictionary of dictionaries, one dictionary per autosome with the geneaccession number as key
    AFregDict = {}
    for chromosome in autosomes:
        AFregDict[chromosome] = {}
    with open("/opt/external_drive/oneill/antifungal/AFimmunity_reg_regions.txt", 'r') as REG:
        for line in REG:
            if line[0] != "#":
                line = line.strip().split()
                chrom,start,stop,name = line[0],int(line[1]),int(line[2]),line[3]
                if chrom in autosomes:
                    AFregDict[chrom][name] = [start,stop,name]
    return(AFregDict)

def make_gene_segment_overlap(dictionaries,snpIDdict,antifungals,AFRegDict,autosomes):
    #from plink.lmiss file generate a list of snps that need to be modified
    with open("151105_CADD_v1.3.tsv", "r") as CADD, open("SNPs_overlap_CADD.ann", "w") as outfile:
        for line in CADD:
            if line[0] != "#":
                line = line.strip().split()
                chrom = "chr" + line[0]
                bp = str(line[1])
                if chrom in autosomes:
                    if bp in dictionaries[chrom]:
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
                        outfile.write("\t".join(line) + "\t" + "".join(dictionaries[chrom][bp]) + "\t" + snpIDdict[chrom][bp] + '\t' + ",".join(candReg) + '\t' + ",".join(candGene) + "\n")

def make_segDict(autosomes):
    #make dictionary of segments
    segDict= {}
    for chromosome in autosomes:
        segDict[chromosome] = {}
    with open("int.overlap", 'r') as CON:
        for line in CON:
            line = line.strip().split()
            seg,ind,startSEG,stopSEG = line[0],str(line[3]),str(line[7]),str(line[8])
            chrom = "chr" + line[4]
            if chrom in autosomes:
                segDict[chrom][seg] = [startSEG,stopSEG,seg,ind]
    return(segDict)

def ann_CADD_seg(segDict,autosomes):
    #from plink.lmiss file generate a list of snps that need to be modified
    with open("SNPs_overlap_CADD.ann", "r") as SNPs, open("final.ann", "w") as outfile:
        for line in SNPs:
            if line[0] != "#":
                line = line.strip().split('\t')
                chrom = "chr" + line[0]
                bp = line[1]
                if chrom in autosomes:
                    for v in segDict[chrom].values():
                        if int(bp) >= int(v[0]) and int(bp) <= int(v[1]):
                            outfile.write("\t".join(line) + "\t" + "\t".join(v) + "\n")



dictionaries,snpIDdict = make_chrdict(autosomes)
print(dictionaries["chr1"]["10250"])
antifungals = make_antifungaldict(autosomes)
print(antifungals["chr8"]["IKBKB"])
AFregDict = make_AFregDict(autosomes)
make_gene_segment_overlap(dictionaries,snpIDdict,antifungals,AFregDict,autosomes)
segDict = make_segDict(autosomes)
for chrom in autosomes:
    print(len(segDict[chrom]))
#!sort -u SNPs_overlap_CADD.ann > SNPs_overlap_CADD_uniq.ann
ann_CADD_seg(segDict,autosomes)
