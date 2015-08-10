setwd('/home/jun/melanoma/')

#####################################################################################

data <- read.csv('Table_S1D.csv', na.strings = c("-", 'n/a', '[Not Available]', '[ERROR]',
                                                 '[Not Applicable]'))
summary(data)
data$Name <- as.character(data$Name)
data$UV.signature
data$MUTATIONSUBTYPES

data$Name[duplicated(gsub('-0[16]{1}$', '', data$Name))]

data <- data[data$Name != "TCGA-ER-A19T-06",]
data <- data[data$Name != "TCGA-ER-A2NF-06",]

ID_BRAF <- data$Name[data$MUTATIONSUBTYPES == 'BRAF_Hotspot_Mutants']
ID_RAS <- data$Name[data$MUTATIONSUBTYPES == 'RAS_Hotspot_Mutants']
ID_NF1 <- data$Name[data$MUTATIONSUBTYPES == 'NF1_Any_Mutants']
ID_TW_UV <- data$Name[data$MUTATIONSUBTYPES == 'Triple_WT' & 
                        data$UV.signature == 'UV signature']
ID_TW_NUV <- data$Name[data$MUTATIONSUBTYPES == 'Triple_WT' & 
                         data$UV.signature == 'not UV']
ID_UV <- data$Name[data$UV.signature == 'UV signature']
ID_NUV <- data$Name[data$UV.signature == 'not UV']

ID_BRAF <- ID_BRAF[!is.na(ID_BRAF)]
ID_RAS <- ID_RAS[!is.na(ID_RAS)]
ID_NF1 <- ID_NF1[!is.na(ID_NF1)]
ID_TW_UV <- ID_TW_UV[!is.na(ID_TW_UV)]
ID_TW_NUV <- ID_TW_NUV[!is.na(ID_TW_NUV)]
ID_UV <- ID_UV[!is.na(ID_UV)]
ID_NUV <- ID_NUV[!is.na(ID_NUV)]


colnames(data)

cnv <- read.delim('./gdac.broadinstitute.org_SKCM.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2013111400.0.0/SKCM.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt')

summary(cnv)
-01D-A192-01
cnv$ID <- as.character(gsub('.-...-....-..$', '' , cnv$Sample))

cnv_BRAF <- cnv[cnv$ID %in% ID_BRAF,1:6]
cnv_RAS <- cnv[cnv$ID %in% ID_RAS,1:6]
cnv_NF1 <- cnv[cnv$ID %in% ID_NF1,1:6]
cnv_TW_UV <- cnv[cnv$ID %in% ID_TW_UV,1:6]
cnv_TW_NUV <- cnv[cnv$ID %in% ID_TW_NUV,1:6]
cnv_UV <- cnv[cnv$ID %in% ID_UV,1:6]
cnv_NUV <- cnv[cnv$ID %in% ID_NUV,1:6]

write.table(cnv_BRAF, file = "CNV_BRAF.txt", append = FALSE, quote = FALSE,
            sep = "\t", dec = ".", row.names = FALSE,
            col.names = TRUE)

write.table(cnv_RAS, file = "CNV_RAS.txt", append = FALSE, quote = FALSE,
            sep = "\t", dec = ".", row.names = FALSE,
            col.names = TRUE)

write.table(cnv_NF1, file = "CNV_NF1.txt", append = FALSE, quote = FALSE,
            sep = "\t", dec = ".", row.names = FALSE,
            col.names = TRUE)

write.table(cnv_TW_UV, file = "CNV_TW_UV.txt", append = FALSE, quote = FALSE,
            sep = "\t", dec = ".", row.names = FALSE,
            col.names = TRUE)

write.table(cnv_TW_NUV, file = "CNV_TW_NUV.txt", append = FALSE, quote = FALSE,
            sep = "\t", dec = ".", row.names = FALSE,
            col.names = TRUE)

write.table(cnv_UV, file = "CNV_UV.txt", append = FALSE, quote = FALSE,
            sep = "\t", dec = ".", row.names = FALSE,
            col.names = TRUE)

write.table(cnv_NUV, file = "CNV_NUV.txt", append = FALSE, quote = FALSE,
            sep = "\t", dec = ".", row.names = FALSE,
            col.names = TRUE)

##########################################################################################

untar('broad.mit.edu_SKCM.Genome_Wide_SNP_6.Level_3.tar.gz', 
      exdir = '.')
con <- file('./broad.mit.edu_SKCM.IlluminaGA_DNASeq.Level_2.1.5.0/skcm_clean_pairs.aggregated.capture.tcga.uuid.somatic.maf')

mut <- read.delim('./broad.mit.edu_SKCM.IlluminaGA_DNASeq.Level_2.1.5.0/skcm_clean_pairs.aggregated.capture.tcga.uuid.somatic.maf',
                  header = FALSE, skip = 5)

colnames(mut) <- unlist(strsplit(readLines(con, n = 5)[5], '\t'))
summary(mut)

head(mut)

pre_seq <- substr(x = mut$ref_context, start = 9, stop = 9)
colnames(mut)
subst <- factor(gsub('^g\\.chr.{0,}:[0-9]{0,}_?[0-9]{0,}','', mut$Genome_Change))
CtoT <- toupper(pre_seq[grep('C>T', subst)]) %in% c('C', 'T')
GtoA <- toupper(pre_seq[grep('G>A', subst)]) %in% c('G', 'A')
CCtoTT <- subst %in% c('CC>TT', 'GG>AA')

sum((CtoT+GtoA+CCtoTT) > 0)


####################################################################################
library('cgdsr')
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)
studies <- getCancerStudies(mycgds)
mycancerstudy = getCancerStudies(mycgds)[78,1]
caselist = getCaseLists(mycgds,mycancerstudy)
mycaselist = getCaseLists(mycgds,mycancerstudy)[2,1]
myclinicaldata = getClinicalData(mycgds,mycaselist)