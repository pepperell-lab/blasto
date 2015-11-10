###combined_filtering_SNPs.py
Script made during the pilot phase to pull various data from a number of files 
and append information to the CADD annotation of variants.

Usage: combined_filtering_SNPs.py

###annotate_eQTLs.py
Takes any file where the first two columns correspond to the chromosome and position
of SNPs (designed for CADD or appended CADD annotation file) and draws info from
a serious of files detailing significant eQTLs from the GTEx portal. Appends to the
input file the tissue and gene if the SNP was found to be an eQTL.

Usage: annotate_eQTLs.py -i [SNPfile] -d [directory of GTEx files]

###ROHsnps_CADD_GTEx.R
A script to filter through output file from annotate_eQTLs. Meant for data exploration,
not as an executable to produce another output.
