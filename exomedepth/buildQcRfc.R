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

}
for (csv in args[2:length(args)]) {
    ed<-rbind(ed,read.table(csv,row.names=1,header=T,fill=TRUE))
}

# remove unclassified rows
ed[ed==""]<-NA
ed<-ed[which(!is.na(ed$status)),]

# factorise status row
ed$status<-as.factor(ed$status)

# build classifier and return error rate
rfc<-randomForest(status~., data=ed, ntree=200, mtry=3, importance=TRUE)

# show result
importance(rfc)
print(rfc)


# define 


# save model
save(list=c("rfc"), file=args[1])
