library(vcfR)

vcf <- vcfR::read.vcfR("C:/Users/Mary/Desktop/sample.vcf")

dat <- vcfR2tidy(vcf, info_only = FALSE, single_frame = FALSE,
          toss_INFO_column = TRUE)

extract_info_tidy(x, info_fields = NULL, info_types = TRUE,
                  info_sep = ";")

extract_gt_tidy(x, format_fields = NULL, format_types = TRUE,
                dot_is_NA = TRUE, alleles = TRUE, allele.sep = "/",
                gt_column_prepend = "gt_", verbose = TRUE)

vcf_field_names(x, tag = "INFO")

#########


##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|GENE_PHENO|SIFT|PolyPhen|DOMAINS|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

dat <- read.table("C:/Users/Mary/Desktop/hmong_roh_regions_VEPannotation_GTEx_antifungal.txt", header=F)
colnames(dat) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", 
                   "BDH001","BDH002","BDH003","BDH004","BDH005","BDH006","BDH007",
                   "BDH008","BDH009","BDH010","BDH018","NA12878",
                   "GTEx","CandReg","CandGene")


#x <- read.table("C:/Users/Mary/Desktop/sample.vcf", skip=4, header=F)
#colnames(x) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

INTO = c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE","STRAND","FLAGS","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID","CANONICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","REFSEQ_MATCH","SOURCE","GENE_PHENO","SIFT","PolyPhen","DOMAINS","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE")
require(tidyr)
require(dplyr)

y <- dat %>% separate(col=INFO, into=INTO, sep = "\\|")
