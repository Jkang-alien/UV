setwd('/home/jun/UV_git/')

#####################################################################################

data <- read.csv('Table_S1D.csv', na.strings = c("-", 'n/a', '[Not Available]', '[ERROR]',
                                                 '[Not Applicable]'))
summary(data)
data$Name <- as.character(data$Name)
data$UV.signature
data$MUTATIONSUBTYPES
data$NECROSIS

data$Name[duplicated(gsub('-0[16]{1}$', '', data$Name))]

data <- data[data$Name != "TCGA-ER-A19T-06",]
data <- data[data$Name != "TCGA-ER-A2NF-06",]

data$chromothripsis [data$SHATTERSEEK_Chromothripsis_calls == 'negative'] <- 'negative'
data$chromothripsis [data$SHATTERSEEK_Chromothripsis_calls != 'negative'] <- 'positive'
data$chromothripsis [is.na(data$SHATTERSEEK_Chromothripsis_calls)] <- NA
data$chromothripsis <- factor(data$chromothripsis)

data$PIGMENT.SCORE.integ <- as.integer(data$PIGMENT.SCORE)
data$LYMPHOCYTE.DISTRIBUTION.integ <- as.integer(data$LYMPHOCYTE.DISTRIBUTION)
data$LYMPHOCYTE.DENSITY.integ <- as.integer(data$LYMPHOCYTE.DENSITY)
data$LYMPHOCYTE.SCORE.integ <- as.integer(data$LYMPHOCYTE.SCORE)


data$PIGMENT.SCORE.f <- factor(data$PIGMENT.SCORE)
data$LYMPHOCYTE.DISTRIBUTION.f <- factor(data$LYMPHOCYTE.DISTRIBUTION)
data$LYMPHOCYTE.DENSITY.f <- factor(data$LYMPHOCYTE.DENSITY)
data$LYMPHOCYTE.SCORE.f <- factor(data$LYMPHOCYTE.SCORE)
data$necrosis <- factor(data$NECROSIS >0 )

########################## TSS #################################]

setwd('/home/jun/melanoma_UV/')
read.delim('tissueSourceSite.txt') -> tss
read.csv('tss.txt') -> location

data$ID <- gsub('-..$', '' , data$Name)
tss[tss$Study.Name == 'Skin Cutaneous Melanoma' , 1:2] -> tss_melanoma
data$tss_ID = as.factor(unlist(strsplit(as.character(data$ID) , '-'))[seq(2,3*dim(data)[1],3)])
data$tss = factor (data$tss_ID, 
                       levels = tss_melanoma[tss_melanoma$TSS.Code %in%  data$tss_ID, ]$TSS.Code,
                       labels = tss_melanoma[tss_melanoma$TSS.Code %in%  data$tss_ID, ]$Source.Site)

data = merge(data, location, by.x = 'tss', by.y = 'Institute', all.x = TRUE)

data$Latitude = round(abs(data$latitude), 0)


##################### Lentigerous #######################

setwd('/home/jun/melanoma_UV')
histologic_type <- read.delim('MELANOMA.txt')
histologic_type$ID <- paste('TCGA-', 
                            gsub('-[A-Z0-9]{3}-[0-9]{2}-[A-Z0-9]{3}', '', histologic_type$ID),
                            sep = '')
summary(histologic_type)

data <- merge(histologic_type, data, by = 'ID', all.y =TRUE)

t(t(colnames(data)))

####################################################################################

library(compareGroups)

compareGroups(UV.signature~ 
                CURATED_AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS +
                CURATED_AGE_AT_TCGA_SPECIMEN +
                GENDER+
                Latitude+
                CURATED_SITE_OF_PRIMARY_TUMOR_KNOWN_PRIMARY_ONLY +
                REGIONAL_VS_PRIMARY+
                type + 
                CURATED_BRESLOW +
                CURATED_ULCERATION +
                PIGMENT.SCORE.integ +
                CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE +
                MUTATIONSUBTYPES+
                NECROSIS+
                necrosis +
                Tumour.content........nuceli.that.are.tumour.cells......0.100..+                         
                LYMPHOCYTE.DISTRIBUTION.integ+                                                                 
                LYMPHOCYTE.DENSITY.integ+                                                                      
                LYMPHOCYTE.SCORE.integ+                                                                      
                PURITY..ABSOLUTE.+                                                                      
                PLOIDY..ABSOLUTE.+                                                                       
                #TOTAL.MUTATIONS+                                                                         
                #UV.RATE+
                chromothripsis ,
              data = data) -> table

table(data$Pathologic_stage, data$uv_g)

createTable(table) -> table

setwd('/home/jun/UV_git/')
export2csv(table, 'table.csv')

compareGroups(  ~ 
                CURATED_AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS +
                CURATED_AGE_AT_TCGA_SPECIMEN +
                GENDER+
                Latitude+
                CURATED_SITE_OF_PRIMARY_TUMOR_KNOWN_PRIMARY_ONLY +
                REGIONAL_VS_PRIMARY+
                type + 
                CURATED_BRESLOW +
                CURATED_ULCERATION +
                PIGMENT.SCORE.integ +
                CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE +
                MUTATIONSUBTYPES+
                NECROSIS+
                necrosis +
                Tumour.content........nuceli.that.are.tumour.cells......0.100..+                         
                LYMPHOCYTE.DISTRIBUTION.integ+                                                                 
                LYMPHOCYTE.DENSITY.integ+                                                                      
                LYMPHOCYTE.SCORE.integ+                                                                      
                PURITY..ABSOLUTE.+                                                                      
                PLOIDY..ABSOLUTE.+                                                                       
                #TOTAL.MUTATIONS+                                                                         
                #UV.RATE+
                chromothripsis ,
              data = data[is.na(data$UV.signature) == FALSE,]) -> table_all

table(data$Pathologic_stage, data$uv_g)

createTable(table_all) -> table_all

setwd('/home/jun/UV_git/')
export2csv(table_all, 'table_all.csv')


chisq.test(data[is.na(data$UV.signature) == FALSE,]$necrosis, data[is.na(data$UV.signature) == FALSE,]$UV.signature)
wilcox.test(CURATED_BRESLOW~UV.signature, data)
wilcox.test(PIGMENT.SCORE.integ~UV.signature, data)
wilcox.test(LYMPHOCYTE.DISTRIBUTION.integ~UV.signature, data)
wilcox.test(LYMPHOCYTE.DENSITY.integ~UV.signature, data)
wilcox.test(LYMPHOCYTE.SCORE.integ~UV.signature, data)
wilcox.test(NECROSIS~UV.signature, data)
wilcox.test(as.numeric(data$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE)~
              data$UV.signature)
wilcox.test(data$Tumour.content........nuceli.that.are.tumour.cells......0.100.. ~ data$UV.signature)


tapply(data$NECROSIS,data$UV.signature, function(x) median(x, na.rm = T) )
tapply(data$Latitude,data$UV.signature, function(x) quantile(x, na.rm = T) )
tapply(data$CURATED_BRESLOW,data$UV.signature, function(x) quantile(x, na.rm = T) )
tapply(data$Tumour.content........nuceli.that.are.tumour.cells......0.100..,data$UV.signature, function(x) quantile(x, na.rm = T) )
tapply(data$NECROSIS,data$UV.signature, function(x) quantile(x, na.rm = T) )

hist(data$LYMPHOCYTE.DISTRIBUTION)
hist(data$LYMPHOCYTE.DENSITY)
hist(data$LYMPHOCYTE.SCORE)

t(t(colnames(data)))

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

CairoSVG(file = "Survival.svg",  width = 5, height = 5, 
         onefile = TRUE, bg = "transparent",
         pointsize = 12)
par(
  #mfrow = c(1,2), 
  mar=c(6,4,2,6), mgp = c(2, 1, 0))


fit = npsurv(Surv(TCGA_M, CURATED_VITAL_STATUS == 'Dead')~
               UV.signature, data = data)
fit

diff = survdiff(Surv(TCGA_M, CURATED_VITAL_STATUS == 'Dead')~
                  UV.signature, data = data)
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
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.3, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(55, 0.87, 'P-value = 0.081', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

dev.off()


################### Multivariate ######################

data$UV.signature <- C(data$UV.signature, contr.treatment, base = 2)

cox_TCGA <- coxph(Surv(TCGA_M, CURATED_VITAL_STATUS == 'Dead')~
                    CURATED_AGE_AT_TCGA_SPECIMEN +
                    #stage (Stage at diagnosis is not compatible with post accession survival)
                    UV.signature, 
                  data = data)

survfit(cox_TCGA)

sink('result.txt')
summary(cox_TCGA)
sink()

############# UV mutation plot ##########################

setwd('/home/jun/UV_git')
data$dipyr_CtoT <- as.numeric(gsub('%', '',data$DIPYRIM.C.T.nTotal.Mut, data$UV.signature))/100

t(t(colnames(data)))

library(Cairo)

CairoSVG(file = "mutation.svg",  width = 5, height = 5, 
         onefile = TRUE, bg = "transparent",
         pointsize = 12)
par(mar=c(4.2,6,4,4))
plot(log(data$TOTAL.MUTATIONS, base = 10), data$dipyr_CtoT^5,
     col =data$UV.signature,
     ylab = 'Proportion of C to T\n at dipyrimidine site',
     xlab = 'Number of mutations',
     xlim = c(1, 5),
     yaxt="n",
     xaxt = 'n'
     )

axis(2, at=c(0.1, 0.6, 0.7, 0.8, 0.9)^5, 
     labels=c(0.1, 0.6, 0.7, 0.8, 0.9), las=2)
axis(1, at=log(c(10, 100, 1000, 10000, 100000), base = 10), 
     labels=c(10, 100, 1000, 10000, 100000), las=1)

legend(1.0, 0.72,
       c('UV signature',
         'not UV signature'),
       pch = 1, 
       col = c(2,1),
       bty = 'n')

dev.off()

library(sm)
# plot densities 
sm.density.compare(data$dipyr_CtoT[is.na(data$dipyr_CtoT) == FALSE]^5, 
                   data$UV.signature[is.na(data$dipyr_CtoT) == FALSE],
                   lty = c(1,1),
                   col = c(1,2),
                   xlim = c(0,0.7),
                   h = 0.05,
                   model = 'none',
                   xlab = '')

sm.density.compare(log(data$TOTAL.MUTATIONS[is.na(data$TOTAL.MUTATIONS) == FALSE],base = 10), 
                   data$UV.signature[is.na(data$TOTAL.MUTATIONS) == FALSE],
                   lty = c(1,1),
                   col = c(1,2),
                   xlim = c(1,5),
                   h = 0.2,
                   model = 'none',
                   xlab = '')
title(main="MPG Distribution by Car Cylinders")

# add legend via mouse click
colfill<-c(2,1) 


data$dipyr_CtoT <- as.numeric(gsub('%', '',data$DIPYRIM.C.T.nTotal.Mut, data$UV.signature))/100

