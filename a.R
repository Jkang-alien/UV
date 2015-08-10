setwd('/home/jun/UV_git/')

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

library(compareGroups)

compareGroups(UV.signature~ 
                CURATED_AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS +
                CURATED_AGE_AT_TCGA_SPECIMEN +
                GENDER+
                CURATED_SITE_OF_PRIMARY_TUMOR_KNOWN_PRIMARY_ONLY +
                REGIONAL_VS_PRIMARY+
                CURATED_BRESLOW +
                CURATED_ULCERATION +
                PIGMENT.SCORE+
                CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE +
                MUTATIONSUBTYPES,
              data = data) -> table

table(data$Pathologic_stage, data$uv_g)

createTable(table) -> table
export2csv(table, 'table.csv')


############################  survival ############################# 
library(survival)
library(rms)

data$stage <- 
  factor(data$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE %in% 
           c('Stage III' ,'Stage IV'),
         levels = c(FALSE, TRUE),
         labels = c('Stage 0-II', 'Stage III-IV'))
data$stage[is.na(data$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE)] <- NA
summary (data$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE)

data$OS_M <- as.numeric(as.character(data$CURATED_DAYS_TO_DEATH_OR_LAST_FU))/30.4
data$TCGA_M <- as.numeric(as.character(data$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU))/30.4

data_early <- subset(data, stage == 'Stage 0-II')
data_late <- subset(data, stage == 'Stage III-IV')

strata = levels(data$UV.signature)

library(Cairo)
CairoSVG(file = "Survival.svg",  width = 10, height = 15, 
         onefile = TRUE, bg = "transparent",
         pointsize = 12)
par(mfrow = c(3,2), mar=c(6,4,2,8), mgp = c(2, 1, 0))

fit = npsurv(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
               UV.signature, data = data)
fit

diff = survdiff(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
                  UV.signature, data = data)
diff

survplot(fit,
         time.inc = 12,
         title = 'All_OA',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.25, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.85, 'P-value = 0.073', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.8, 'All OA', cex = 1,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

fit = npsurv(Surv(TCGA_M, CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0....ALIVE.OR.CENSORED...1....DEAD.OF.MELANOMA.. == 1)
             ~ UV.signature, data = data)
fit

diff = survdiff(Surv(TCGA_M, CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0....ALIVE.OR.CENSORED...1....DEAD.OF.MELANOMA.. == 1)
                ~ UV.signature, data = data)
diff

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.25, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.85, 'P-value = 0.33', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.8, 'All TCGA', cex = 1,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

################## Early Stage ##############################
fit = npsurv(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
               UV.signature, data = data_early)
fit

diff = survdiff(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
                  UV.signature, data = data_early)
diff

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.25, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.85, 'P-value = 0.486', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.8, 'Early stage Overall', cex = 1,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

fit = npsurv(Surv(TCGA_M, CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0....ALIVE.OR.CENSORED...1....DEAD.OF.MELANOMA.. == 1)
             ~ UV.signature, data = data_early)
fit

diff = survdiff(Surv(TCGA_M, CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0....ALIVE.OR.CENSORED...1....DEAD.OF.MELANOMA.. == 1)
                ~ UV.signature, data = data_early)
diff

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.25, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.85, 'P-value = 0.645', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.8, 'Early stage TCGA', cex = 1,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

###################### Late Stage ########################

fit = npsurv(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
               UV.signature, data = data_late)
fit

diff = survdiff(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
                  UV.signature, data = data_late)
diff

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.25, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.85, 'P-value = 0.046', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.8, 'Late stage Overall', cex = 1,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

fit = npsurv(Surv(TCGA_M, CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0....ALIVE.OR.CENSORED...1....DEAD.OF.MELANOMA.. == 1)
             ~ UV.signature, data = data_late)
fit

diff = survdiff(Surv(TCGA_M, CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0....ALIVE.OR.CENSORED...1....DEAD.OF.MELANOMA.. == 1)
                ~ UV.signature, data = data_late)
diff

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:3),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.25, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.85, 'P-value = 0.099', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.8, 'Late stage TCGA', cex = 1,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')
dev.off()
