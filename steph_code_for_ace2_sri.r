#code from Steph she used for ACE2 Isoform analysis

library(ASpli)
library(GenomicFeatures)
setwd("C:/Users/snerella/Box Sync/ucsf/Spiromics")
TxDb <- makeTxDbFromGFF(
  file="il6r_fas.gtf",
  format="gtf")

files<-list.files("il6r_fas_bam_subsection/")
                  #C:/Users/christensons/Dropbox/Lab/Projects/Spiromics/EpithelialSequencing/SPIROMICS ACE2 Analysis/ACE2_bam_subsection/ace2_bam_subsection/SAMtoBAM")

bamfiles<-files[seq(1, length(files), 2)]
bamfiles<-bamfiles[-104]
indexfiles<-files[seq(2, length(files), 2)]
indexfiles<-indexfiles[-104]
#Sample and Plate tracking for 2017 Run
p1<-read.csv("SPIROMICS_EpiSeq_Plate1Tracking.csv", fileEncoding = "latin1")
p2<-read.csv("SPIROMICS_EpiSeq_Plate2Tracking.csv", fileEncoding = "latin1")
colnames(p2)<-colnames(p1)

p1$plate<-1
p2$plate<-2
pt<-rbind(p1, p2)

pt$SAMID<-sapply(as.character(pt$Full.RNA.Sample.ID),function(x) strsplit(x,"_644__")[[1]][2])
pt$SAMID<-sapply(as.character(pt$SAMID),function(x) strsplit(x,"__")[[1]][1])

samples<-as.list(bamfiles)
sampledata<-as.matrix(bamfiles)
samples = sapply(as.character(samples),function(x) strsplit(x,"_646__")[[1]][2])
samples[1:106]<-names(samples[1:106])
samples = sapply(as.character(samples),function(x) strsplit(x,"_il6r")[[1]][1])
samples[83:106]<-sapply(as.character(samples[83:106]),function(x) strsplit(x,"_")[[1]][3])

samplenames<-as.data.frame(cbind(bamfiles, samples))
colnames(samplenames)<-c("colnames", "samid")

Run1_IDmatch <- read.csv("GenentechSPIROMICSEpiSeqRun_SampleIDMatch.csv",row.names=1)
colnames(Run1_IDmatch)<-c('samid', 'bronchid', 'subjid', 'plate')
R1SubjIDs <- as.character(Run1_IDmatch$subjid[match(samplenames$samid, Run1_IDmatch$samid)])
R1SubjIDs[83:106]<-samples[83:106]
plate <- as.numeric(Run1_IDmatch$plate[match(samplenames$samid, Run1_IDmatch$samid)])
plate[83:106]<-3
samplenames$subjid<-R1SubjIDs
samplenames$plate<-plate


samplenames$subjid[duplicated(samplenames$subjid)]
grep("WF122906",samplenames$subjid)#remove 65
grep("UT172747",samplenames$subjid)#remove 143
grep("CU104873",samplenames$subjid)#remove 1

samplenames<-samplenames[-143,]
samplenames<-samplenames[-65,]
samplenames<-samplenames[-1,]
core6metadata <- read.csv("C:/Users/christensons/Dropbox/Lab/Projects/Spiromics/Core6/CORE6_1_clinical_191118.csv")
bronchmeta<-read.csv("C:/Users/christensons/Dropbox/Lab/Projects/Spiromics/EpithelialSequencing/SPIROMICS ACE2 Analysis/SPIROMICSbronchmetadata.csv")
samplenames$smoking<-as.factor(core6metadata$CURRENT_SMOKER_V1[match(samplenames$subjid, core6metadata$SUBJID)])
samplenames$age<-as.numeric(core6metadata$AGE_DERV_V1[match(samplenames$subjid, core6metadata$SUBJID)])
samplenames$gender<-as.factor(core6metadata$GENDER[match(samplenames$subjid, core6metadata$SUBJID)])
samplenames$BMI<-as.numeric(core6metadata$BMI_CM_V1[match(samplenames$subjid, core6metadata$SUBJID)])
samplenames$race<-as.factor(core6metadata$RACE[match(samplenames$subjid, core6metadata$SUBJID)])
samplenames$obesity<-NA
samplenames$obesity[samplenames$BMI>29.99]<-"yes"
samplenames$obesity[samplenames$BMI<30]<-"no"
samplenames$obesity<-as.factor(samplenames$obesity)
samplenames$cvdisease<-as.factor(core6metadata$CARDIOVASCULAR_CONDITION_V1_A[match(samplenames$subjid, core6metadata$SUBJID)])
samplenames$HTN<-as.factor(core6metadata$HYPERTENSION_DERV_V1[match(samplenames$subjid, core6metadata$SUBJID)])
samplenames$sleepapnea<-as.factor(core6metadata$SLEEP_APNEA_DERV_V1[match(samplenames$subjid, core6metadata$SUBJID)])
samplenames$diabetes<-as.factor(core6metadata$DIABETES_DERV_V1[match(samplenames$subjid, core6metadata$SUBJID)])
samplenames$ifn<-as.numeric(bronchmeta$ifncenteredmean[match(samplenames$subjid, bronchmeta$SUBJID)])
median(samplenames$ifn, na.rm=TRUE)
samplenames$ifncat<-NA
samplenames$ifncat[samplenames$ifn<0.001]<-"low"
samplenames$ifncat[samplenames$ifn>0]<-"high"
samplenames$ifncat<-as.factor(samplenames$ifncat)

bamfiles2<-bamfiles[match(samplenames$colnames, bamfiles)]
targets <- data.frame( row.names = samplenames$subjid,
                       bam = bamfiles2,
                       smoking = samplenames$smoking,
                       stringsAsFactors = FALSE )

# extract features from annotation
features <- binGenome( TxDb )
# Accesors of ASpliFeatures object
geneCoord <- featuresg( features )
binCoord <- featuresb( features )
junctionCoord <- featuresj( features )
# Extract metadata annotation, such as bins names, event and feature type
binMetadata <- mcols( binCoord )
setwd("C:/Users/christensons/Dropbox/Lab/Projects/Spiromics/EpithelialSequencing/SPIROMICS ACE2 Analysis/ACE2_bam_subsection/ace2_bam_subsection/SAMtoBAM")
bam <- loadBAM(targets)
counts <- readCounts (
  features,
  bam,
  targets,
  cores = 1,
  readLength = 150L,
  maxISize = 50000,
  minAnchor = 10 )
# Accessing count and read density data
GeneCounts <- countsg(counts)
GeneRd <- rdsg(counts)
BinCounts <- countsb(counts)
BinRd <- rdsb(counts)
JunctionCounts <- countsj(counts)
# Export count data to text files
write.csv(BinCounts, "ACE2_Exon1_Counts.csv")
write.csv(BinRd, "ACE2_Exon1_Reads.csv")
targets$smoking[54]<-0
targets$smoking[87]<-0
du <- DUreportBinSplice( counts,
                         targets,
                         minGenReads = 0,
                         minBinReads = 0,
                         minRds = 0,
                         contrast = NULL,
                         forceGLM = FALSE,
                         ignoreExternal = TRUE,
                         ignoreIo = TRUE,
                         ignoreI = TRUE,
                         filterWithContrasted = FALSE)
#bam <- scanBam("C:/Users/christensons/Dropbox/Lab/Projects/Spiromics/EpithelialSequencing/SPIROMICS ACE2 Analysis/ACE2_bam_subsection/ace2_bam_subsection/SAMtoBAM/SAM24314828_ace2_region.bam.bam")


ExonCts<-BinCounts[,-(1:9)]
ExonCtsT<-as.data.frame(t(ExonCts))
ExonCtsT<-ExonCtsT[,-2]
ExonCtsT$dACE2<-ExonCtsT$`ENSG00000130234:E001`/.140
ExonCtsT$Exon1b<-ExonCtsT$`ENSG00000130234:E002`/.289
ExonCtsT$Exon1a<-ExonCtsT$`ENSG00000130234:E003`/.203
#ExonCtsT$sum<-ExonCtsT$`ENSG00000130234:E001`+ ExonCtsT$`ENSG00000130234:E002`+ 
ExonCtsT$`ENSG00000130234:E003`

sfs<-read.csv("SizeFactorsForExonAnalysis.csv", row.names=1)
ExonCtsT$sizefactors<-sfs$x[match(row.names(ExonCtsT), row.names(sfs))]
ExonCtsTmin<-ExonCtsT[-which(is.na(ExonCtsT$sizefactors)),]
ExonCtsTmin$dACE2norm<-ExonCtsTmin$dACE2*ExonCtsTmin$sizefactors
ExonCtsTmin$Exon1bnorm<-ExonCtsTmin$Exon1b*ExonCtsTmin$sizefactors
ExonCtsTmin$Exon1anorm<-ExonCtsTmin$Exon1a*ExonCtsTmin$sizefactors
ExonCtsTmin$Exon1abnorm<-(ExonCtsTmin$Exon1a+ ExonCtsTmin$Exon1b)*ExonCtsTmin$sizefactors

NormExon<-ExonCtsTmin[,(8:11)]
NormExon[NormExon==0]<-1
ExonCtsTmin[,(8:11)]<-NormExon
ExonCtsTmin$dACE2log2norm<-log2(ExonCtsTmin$dACE2norm)
ExonCtsTmin$Exon1blog2norm<-log2(ExonCtsTmin$Exon1bnorm)
ExonCtsTmin$Exon1alog2norm<-log2(ExonCtsTmin$Exon1anorm)
ExonCtsTmin$Exon1ablog2norm<-log2(ExonCtsTmin$Exon1abnorm)
write.csv(ExonCtsTmin, "Log2NormalizedACE2Exon1.csv") 

library(psych)
ExonCtsTmin$geomean<-apply<- apply(ExonCtsTmin[,(1:3)], 1, geometric.mean)
ExonCtsTmin$dACE2geonorm<-ExonCtsTmin$dACE2*ExonCtsTmin$geomean
ExonCtsTmin$Exon1bgeonorm<-ExonCtsTmin$Exon1b*ExonCtsTmin$geomean
ExonCtsTmin$Exon1ageonorm<-ExonCtsTmin$Exon1a*ExonCtsTmin$geomean
ExonCtsTmin$Exon1abgeonorm<-(ExonCtsTmin$Exon1a+ExonCtsTmin$Exon1b)*ExonCtsTmin$geomean
#ExonCtsT<-ExonCtsT[-153,]
NormExon<-ExonCtsTmin[,(17:20)]
NormExon[NormExon==0]<-1
ExonCtsTmin[,(17:20)]<-NormExon
ExonCtsTmin$dACE2log2geonorm<-log2(ExonCtsTmin$dACE2geonorm)
ExonCtsTmin$Exon1blog2geonorm<-log2(ExonCtsTmin$Exon1bgeonorm)
ExonCtsTmin$Exon1alog2geonorm<-log2(ExonCtsTmin$Exon1ageonorm)
ExonCtsTmin$Exon1ablog2geonorm<-log2(ExonCtsTmin$Exon1abgeonorm)

write.csv(ExonCtsTmin, "Log2NormalizedACE2Exon1.csv") 

ExonCtsTmin$dACE2raw<-ExonCtsTmin$`ENSG00000130234:E001`
ExonCtsTmin$Exon1braw<-ExonCtsTmin$`ENSG00000130234:E002`
ExonCtsTmin$Exon1araw<-ExonCtsTmin$`ENSG00000130234:E003`
ExonCtsTmin$dACE2rawnorm<-ExonCtsTmin$dACE2raw*ExonCtsTmin$sizefactors
ExonCtsTmin$Exon1brawnorm<-ExonCtsTmin$Exon1braw*ExonCtsTmin$sizefactors
ExonCtsTmin$Exon1arawnorm<-ExonCtsTmin$Exon1araw*ExonCtsTmin$sizefactors
ExonCtsTmin$Exon1abrawnorm<-(ExonCtsTmin$Exon1araw + ExonCtsTmin$Exon1braw)*ExonCtsTmin$sizefactors

NormExon<-ExonCtsTmin[,(28:31)]
NormExon[NormExon==0]<-1
ExonCtsTmin[,(28:31)]<-NormExon
ExonCtsTmin$dACE2log2rawnorm<-log2(ExonCtsTmin$dACE2rawnorm)
ExonCtsTmin$Exon1blog2rawnorm<-log2(ExonCtsTmin$Exon1brawnorm)
ExonCtsTmin$Exon1alog2rawnorm<-log2(ExonCtsTmin$Exon1arawnorm)
ExonCtsTmin$Exon1ablog2rawnorm<-log2(ExonCtsTmin$Exon1abrawnorm)

load("C:/Users/christensons/Dropbox/Lab/Projects/RNAseq/SARP_RITA_SPIROMICS_RawCounts/SPIROMICS_compiled/SPIROMICSepi_raw_protcoding_Jan2020.RData")
grep("ENSG00000130234", rownames(SPIROMICS_fullcounts_pc_QA))
ACE2raw<-SPIROMICS_fullcounts_pc_QA[5272,]
ExonCtsTmin$ACE2raw<-ACE2raw[match(rownames(ExonCtsTmin), names(ACE2raw))]
ExonCtsTmin$ACE2rawmean<-mean(ExonCtsTmin$ACE2)
ExonCtsTmin$dACE2genecorr<-((ExonCtsTmin$dACE2raw*ExonCtsTmin$ACE2rawmean)/ExonCtsTmin$ACE2raw)*ExonCtsTmin$sizefactors
ExonCtsTmin$Exon1bgenecorr<-((ExonCtsTmin$Exon1braw*ExonCtsTmin$ACE2rawmean)/ExonCtsTmin$ACE2raw)*ExonCtsTmin$sizefactors
ExonCtsTmin$Exon1agenecorr<-((ExonCtsTmin$Exon1araw*ExonCtsTmin$ACE2rawmean)/ExonCtsTmin$ACE2raw)*ExonCtsTmin$sizefactors
ExonCtsTmin$Exon1abgenecorr<-(((ExonCtsTmin$Exon1araw+ExonCtsTmin$Exon1b)*ExonCtsTmin$ACE2rawmean)/ExonCtsTmin$ACE2raw)*ExonCtsTmin$sizefactors
#ExonCtsT<-ExonCtsT[-153,]
NormExon<-ExonCtsTmin[,(38:41)]
NormExon[NormExon==0]<-1
ExonCtsTmin[,(38:41)]<-NormExon
ExonCtsTmin$dACE2log2genecorr<-log2(ExonCtsTmin$dACE2genecorr)
ExonCtsTmin$Exon1blog2genecorr<-log2(ExonCtsTmin$Exon1bgenecorr)
ExonCtsTmin$Exon1alog2genecorr<-log2(ExonCtsTmin$Exon1agenecorr)
ExonCtsTmin$Exon1ablog2genecorr<-log2(ExonCtsTmin$Exon1abgenecorr)

load("C:/Users/christensons/Dropbox/Lab/Projects/RNAseq/SARP_RITA_SPIROMICS_RawCounts/SPIROMICS_compiled/SPIROMICSepi_rlog_protcoding_Jan2020.RData")
dim(rld_QA_assay)
epidata<-read.csv("C:/Users/christensons/Dropbox/Lab/Projects/Spiromics/EpithelialSequencing/SPIROMICS ACE2 Analysis/ACE2_bam_subsection/ace2_bam_subsection/ACE2IsoformAnalysis/Epibrush_Metadata.csv")
epidata$plate<-as.factor(epidata$plate)

rld_ACE<- rld_QA_assay[,match(rownames(ExonCtsTmin), colnames(rld_QA_assay))]
epidata_ACE<- epidata[match(rownames(ExonCtsTmin), epidata$SUBJID),]
rm(epidata, rld_QA_assay)

rld_ACE2<-rbind(rld_ACE, ExonCtsTmin$dACE2log2norm, ExonCtsTmin$Exon1blog2norm,
                ExonCtsTmin$Exon1alog2norm,ExonCtsTmin$dACE2log2geonorm,
                ExonCtsTmin$Exon1blog2geonorm,ExonCtsTmin$Exon1alog2geonorm,
                ExonCtsTmin$dACE2log2rawnorm,ExonCtsTmin$Exon1blog2rawnorm,
                ExonCtsTmin$Exon1alog2rawnorm,ExonCtsTmin$Exon1ablog2norm,
                ExonCtsTmin$Exon1ablog2geonorm,ExonCtsTmin$Exon1ablog2rawnorm,
                ExonCtsTmin$dACE2log2genecorr,ExonCtsTmin$Exon1blog2genecorr,
                ExonCtsTmin$Exon1alog2genecorr,ExonCtsTmin$Exon1ablog2genecorr)
rownames(rld_ACE2)[14844:14859]<-c('dACE2log2norm','Exon1blog2norm',
                                   'Exon1alog2norm','dACE2log2geonorm',
                                   'Exon1blog2geonorm','Exon1alog2geonorm',
                                   'dACE2log2rawnorm','Exon1blog2rawnorm',
                                   'Exon1alog2rawnorm','Exon1ablog2norm',
                                   'Exon1ablog2geonorm','Exon1ablog2rawnorm',
                                   'dACE2log2genecorr','Exon1blog2genecorr',
                                   'Exon1alog2genecorr','Exon1ablog2genecorr')
library(sva)
plate <- epidata_ACE$plate
modcombat = model.matrix(~1, data=epidata_ACE) #No other adjustment variables
combat_rld_QA <- ComBat(dat=rld_ACE2, batch=plate, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

grep("ENSG00000130234", rownames(combat_rld_QA))
epidata_ACE<-cbind(epidata_ACE, ExonCtsTmin$dACE2log2norm, ExonCtsTmin$Exon1blog2norm,
                   ExonCtsTmin$Exon1alog2norm,ExonCtsTmin$dACE2log2geonorm,
                   ExonCtsTmin$Exon1blog2geonorm,ExonCtsTmin$Exon1alog2geonorm,
                   ExonCtsTmin$dACE2log2rawnorm,ExonCtsTmin$Exon1blog2rawnorm,
                   ExonCtsTmin$Exon1alog2rawnorm,ExonCtsTmin$Exon1ablog2norm,
                   ExonCtsTmin$Exon1ablog2geonorm,ExonCtsTmin$Exon1ablog2rawnorm,
                   ExonCtsTmin$dACE2log2genecorr,ExonCtsTmin$Exon1blog2genecorr,
                   ExonCtsTmin$Exon1alog2genecorr,ExonCtsTmin$Exon1ablog2genecorr)
colnames(epidata_ACE)[708:723]<-c('dACE2log2norm','Exon1blog2norm',
                                  'Exon1alog2norm','dACE2log2geonorm',
                                  'Exon1blog2geonorm','Exon1alog2geonorm',
                                  'dACE2log2rawnorm','Exon1blog2rawnorm',
                                  'Exon1alog2rawnorm','Exon1ablog2norm',
                                  'Exon1ablog2geonorm','Exon1ablog2rawnorm',
                                  'dACE2log2genecorr','Exon1blog2genecorr',
                                  'Exon1alog2genecorr','Exon1ablog2genecorr')
epidata_ACE<-cbind(epidata_ACE, t(combat_rld_QA[14844:14859,]))
colnames(epidata_ACE)[724:739]<-c('dACE2log2normcombat','Exon1blog2normcombat',
                                  'Exon1alog2normcombat','dACE2log2geonormcombat',
                                  'Exon1blog2geonormcombat','Exon1alog2geonormcombat',
                                  'dACE2log2rawnormcombat','Exon1blog2rawnormcombat',
                                  'Exon1alog2rawnormcombat','Exon1ablog2normcombat',
                                  'Exon1ablog2geonormcombat','Exon1ablog2rawnormcombat',
                                  'dACE2log2genecorrcombat','Exon1blog2genecorrcombat',
                                  'Exon1alog2genecorrcombat','Exon1ablog2genecorrcombat')
write.csv(ExonCtsTmin, "ACE2_Exon1_Counts_BeforeBatchCorrection.csv")
write.csv(epidata_ACE, "epidatabronch_withACE2Exonreads_beforeANDafterCOMBAT.csv")
