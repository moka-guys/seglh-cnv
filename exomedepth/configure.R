# configure.R
# -----------------------------------------------------------
# Builds a random forest classifier from QC data (csv output)
# Make sure to amend a column labeled status with known classifications (eg. PASS/FAIL)

# Usage
# -----
# RScript buildQcRfc.R output.RData <csv1> (<csv2> ...)

## install packages
options(menu.graphics=FALSE)
options(repos=structure(c(CRAN="http://cran.ma.imperial.ac.uk/")))
if (!("randomForest" %in% installed.packages())) install.packages("randomForest",dependencies=TRUE)
# biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c('Biostrings','GenomicRanges'))

## load requirements
require(randomForest)
require(GenomicRanges)

## parameters for random forest
ntree<-200
mtrytune<-TRUE

# read command line args
args<-commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop('No output file specified')
}

# --- Build Random Forest Classifier ---
rfc<-NULL
print(args)
csvs<-args[which(endsWith(args,'.csv'))]
beds<-args[which(endsWith(args,'.bed'))]
if (length(csvs)>0) {
  # read data
  ed<-data.frame()
  for (csv in csvs) {
      ed<-rbind(ed,read.table(csv,row.names=1,header=T,fill=TRUE))
  }
  # remove unclassified rows
  ed[ed==""]<-NA
  ed<-ed[which(!is.na(ed$status)),]
  # factorise status row
  ed$status<-as.factor(ed$status)
  # build classifier and return error rate
  variables<-ncol(ed)-1
  oob.err<-double(variables)
  best.mtry<-floor(sqrt(variables))
  if (mtrytune) {
    for (mtry in 1:variables) {
      rf<-randomForest(status~., data=ed, ntree=ntree, mtry=mtry, importance=TRUE)
      oob.err[mtry]=rf$err.rate[nrow(rf$err.rate),"OOB"]
    }
    # pick best classifier
    best.mtry<-which(min(oob.err)==oob.err)
  }
  print(paste("MTRY",best.mtry))
  rfc<-randomForest(status~., data=ed, ntree=ntree, mtry=best.mtry, importance=TRUE)
  # show result
  importance(rfc)
  print(rfc)
}

# --- get annotation bed ofiles ---
annotations<-NA
if (length(beds)>0) {
  message('Getting Segment annotations...')
  ranges<-list()
  for (bed in beds) {
    ranges[[bed]]<-with(read.table(bed,header=FALSE), GRanges(seqnames = V1, IRanges(start=V2+1, end=V3, names=V4),names=gsub("_"," ",V4)))  # BED FILE 0-based
  }
  annotations<-unlist(as(ranges, "GRangesList"))
}
# --- QC threshold overrides (edit as needed) ---
limits<-list(
  medcor=c(NA, 0.90),    # median correlation within batch
  maxcor=c(0.95, 0.90),  # max correlation within batch
  refcor=c(0.95, 0.90),  # reference set correlation
  refcount=c(3,1),       # reference set size (selected reference samples)
  coeffvar=c(120, 150),  # coefficient of variation
  coverage=c(100),       # exon coverage limit
  expectedbf=c(5.0, NA), # expected BF
  minrefs=c(2,Inf)       # minimum reference set size
)

# save model
save(list=c("rfc","limits","annotations"), file=args[1])
