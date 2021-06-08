# R script for exome depth data (takes a readcount data object created by readCount.R)
#
# INSTALL PACKAGES
#
options(menu.graphics=FALSE)
options(repos=structure(c(CRAN="http://cran.ma.imperial.ac.uk/")))
if (!("ExomeDepth" %in% installed.packages())) {
  source("http://bioconductor.org/biocLite.R")
  biocLite('Biostrings',ask=FALSE)
  biocLite('Rsamtools',ask=FALSE)
  biocLite('GenomicRanges',ask=FALSE)
  biocLite('GenomicAlignments',ask=FALSE)
  install.packages("ExomeDepth",dependencies=TRUE)
}

#
# LOAD LIBRARIES
#
cat(paste('Loading Packages\n'))
require(warn.conflicts=FALSE,quietly=FALSE,package="GenomicRanges")
require(warn.conflicts=FALSE,quietly=FALSE,package="ExomeDepth")
require(warn.conflicts=FALSE,quietly=FALSE,package="randomForest")
require(warn.conflicts=FALSE,quietly=FALSE,package="Rsamtools")

#
# READ ARGS
#
# inferNormals.R
#   readCount.RData

args<-commandArgs(trailingOnly = TRUE)
setwd(dirname(args[1]))  # set working directory to target directory
load(args[1])
tries<-ifelse(is.na(args[2]),1,as.numeric(args[2]))
from<-ifelse(is.na(args[3]),1,as.numeric(args[3]))
to<-ifelse(is.na(args[4]),length(testsamplenames),as.numeric(args[4]))

#
# print what will be printed
#

message(paste('Read Count Data:',args[1]))

#
# Check if run in batch mode
#
if (!all(refsamplenames==testsamplenames)) {
  stop(paste("You must prepare data in batch normalisation mode."))
}

#
# run if enough samples
#
results<-list()
cnvcount<-vector()
refcount<-vector()
passed<-vector()
i<-0
bed.data<-data.frame()
ref.data<-data.frame(
  test=c(),
  ref=c()
)

exons<-GRanges(seqnames = rois$chromosome,
  IRanges(start=rois$start,end=rois$end),
  names = rois$name)

score.dups<-vector()
names.dups<-vector()
score.dels<-vector()
names.dels<-vector()
for (testsample in testsamplenames[from:to]) {
  print(paste(i+from,to,sep='/'))
  i<-i+1
  samplecnvs<-data.frame()
  for (t in 1:tries) {
    #
    # pick reference
    #
    refsamples<-refsamplenames[which(refsamplenames!=testsample)]
    if (tries>1) refsamples<-sample(refsamples,floor(length(refsamples)/tries))
    refset<-suppressWarnings(select.reference.set(
      test.counts=counts[,testsample],
      reference.counts=as.matrix(counts[,refsamples]),
      bin.length=(counts$end - counts$start)
    ))
    #
    # prepare reference (sum reference choice)
    #
    reference.selected<-apply(
      X = as.matrix(counts[, refset$reference.choice, drop = FALSE]),
      MAR = 1, FUN = sum)
    reference.selected<-as.vector(reference.selected)
    #
    # run exome depth
    #
    message('*** Creating ExomeDepth object...')
    suppressWarnings(ED <- new('ExomeDepth',
                              test = counts[,testsample],
                              reference = reference.selected,
                              formula = 'cbind(test, reference) ~ 1'))
    #
    # call CNV
    #
    message('*** Calling CNVs...')
    result<-CallCNVs(x = ED,
                    transition.probability = 10^-4,
                    chromosome = counts$chromosome,
                    start = counts$start, 
                    end =counts$end,
                    name = counts$exon)
    if (length(result@CNV.calls)>0) {
      cnvs<-result@CNV.calls[result@CNV.calls[,9]>0,]
      if (nrow(cnvs)>0) samplecnvs<-rbind(samplecnvs,cnvs[,c(1,2,9,10,11,12)])
    }
  }
  dups<-table(unlist(apply(samplecnvs[which(samplecnvs[,6]>1),], 1, function(x) x[1]:x[2])))==tries
  score.dups<-append(score.dups,length(names(dups)[which(dups)]))
  names.dups<-append(names.dups,paste(names(dups)[which(dups)],collapse=","))
  dels<-table(unlist(apply(samplecnvs[which(samplecnvs[,6]<1),], 1, function(x) x[1]:x[2])))==tries
  score.dels<-append(score.dels,length(names(dels)[which(dels)]))
  names.dels<-append(names.dels,paste(names(dels)[which(dels)],collapse=","))

  # annotate
  # if (length(result@CNV.calls)>0) {
  #   # add exon numbers (from subset)
  #   result<-AnnotateExtra(
  #     x = result,
  #     reference.annotation = exons,
  #     min.overlap = 0.0001,
  #     column.name = 'exons.hg19')
  #   result@annotations$name<-as.factor(sapply(strsplit(as.character(result@annotations$name),'_'),"[[",1))
  #   # write BED file of CNVs
  #   cnv<-result@CNV.calls[which(result@CNV.calls[,9]>0),]
  #   bed.data<-rbind(bed.data,cbind(cnv,sample=testsample))
  # }
  # filter by sex match
  # sex<-unlist(strsplit(testsample,'_'))[3]
  # refsex<-sapply(strsplit(refsets[[testsample]]$reference.choice,'_'),function(x) x[3])
  # if (all(refsex==sex)) {
  #   results[[testsample]]<-result
  #   cnvcount<-append(cnvcount,nrow(result@CNV.calls[which(result@CNV.calls$BF>10),]))
  #   refcount<-append(refcount,length(refsets[[testsample]]$reference.choice))
  # }
  # for (r in refsets[[testsample]]$reference.choice) {
  #   ref.data<-rbind(ref.data,c(testsample,r))
  # }
}

# stable CNV score
result<-data.frame(
  sample=testsamplenames[from:to],
  score.dups=score.dups,
  score.dels=score.dels,
  names.dups=names.dups,
  names.dels=names.dels)
cnv.file<-sub("[.][^.]*$", ".scores", args[1], perl=TRUE)
write.table(result, file=cnv.file, sep='\t', col.names=FALSE, quote=FALSE, row.names=FALSE, append=TRUE)
# CNV
# write.table(bed.data, file=sub("[.][^.]*$", ".bed", args[1], perl=TRUE), sep='\t', col.names=FALSE, quote=FALSE, row.names=FALSE)
# LINKS
# ref.file<-sub("[.][^.]*$", ".refs", args[1], perl=TRUE)
# write.table(ref.data, file=ref.file, sep='\t', col.names=FALSE, quote=FALSE, row.names=FALSE)
# save image
# data.file<-sub("[.][^.]*$", ".out.RData", args[1], perl=TRUE)
# save.image(data.file)
