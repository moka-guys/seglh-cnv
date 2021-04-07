# buildQcRfc
# -----------------------------------------------------------
# Builds a random forest classifier from QC data (csv output)
# Make sure to amend a column labeled status with known classifications (eg. PASS/FAIL)

# Usage
# -----
# RScript buildQcRfc.R output.RData <csv1> (<csv2> ...)

require(randomForest)

# read command line args
args<-commandArgs(trailingOnly=TRUE)

# --- Build Random Forest Classifier ---
rfc<-NULL
if (length(args)>1) {
  # read data
  ed<-data.frame()
  for (csv in args[2:length(args)]) {
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
  for (mtry in 1:variables) {
    rf<-randomForest(status~., data=ed, ntree=200, mtry=mtry, importance=TRUE)
    oob.err[mtry]=rf$err.rate[nrow(rf$err.rate),"OOB"]
  }
  # pick best classifier
  best.mtry<-which(min(oob.err)==oob.err)
  print(paste("MTRY",best.mtry))
  rfc<-randomForest(status~., data=ed, ntree=200, mtry=best.mtry, importance=TRUE)
  # show result
  importance(rfc)
  print(rfc)
}

# --- QC threshold overrides (edit as needed) ---
limits<-list(
  medcor=c(NA, 0.90),   # median correlation within batch
  maxcor=c(0.95, 0.90), # max correlation within batch
  refcor=c(0.95, 0.90), # reference set correlation
  refcount=c(5,3),      # refernce set size (selected reference samples)
  coeffvar=c(30, 35),   # coefficient of variation
  coverage=c(100)       # exon coverage limit
)

# save model
save(list=c("rfc","limits"), file=args[1])
