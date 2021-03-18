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
if (!("knitr" %in% installed.packages())) {
  install.packages("knitr",dependencies=TRUE)
}

#
# LOAD LIBRARIES
#
cat(paste('Loading Packages\n'))
require(warn.conflicts=FALSE,quietly=FALSE,package="GenomicRanges")
require(warn.conflicts=FALSE,quietly=FALSE,package="ExomeDepth")
require(warn.conflicts=FALSE,quietly=FALSE,package="knitr")

#
# DEFINE FUNCTIONS
#
selectReferenceSet<-function(df, test.sample, ref.samples, remove.test=TRUE) {
    # remove self
    if (remove.test) ref.samples<-ref.samples[which(ref.samples!=test.sample)]
    # select best reference samples
    select.reference.set(test.counts = df[,test.sample], reference.counts = as.matrix(df[,ref.samples]), bin.length = (df$end - df$start))
}

#
# READ ARGS
#
cmd<-commandArgs(trailingOnly = FALSE)
runningScript<-unlist(strsplit(cmd[which(substr(cmd,1,7)=='--file=')],'='))[2]
scriptDirectory<-normalizePath(dirname(runningScript))
args<-commandArgs(trailingOnly = TRUE)
setwd(dirname(args[2]))  # set working directory to target directory (THIS IS ABSOLUTELY NECESSARY!)

#
# print what will be printed
#
panel<-unlist(strsplit(args[3],':'))  # covered regons name
testsample<-basename(unlist(strsplit(args[5],':'))[1])  # RGID
samplename<-as.character(unlist(strsplit(args[5],':'))[2]) # RGSM
threshold<-as.numeric(unlist(strsplit(args[5],':'))[3]) # Threshold set in ngs_config
if(is.na(threshold)) {
    threshold<-10^-4
}
pipeversion<-args[1]
load(args[4])  # loads counts,samplenames,rois (extended exons.hg19)
message(paste('        Version:',args[1]))
message(paste('         Output:',args[2]))
message(paste('            ROI:',panel[1]))
message(paste('     Panel Name:',panel[2]))
message(paste('    Read Counts:',args[4]))
message(paste('   Original BAM:',testsample))
message(paste('     SampleName:',samplename))
message(paste('    Ref samples:',length(bam)-1))

#
# run if enough samples
#
refsets<-list()
results<-list()
if (length(bam)>=3) {
  #
  # read exons/ROI and create subset
  #
  message('Getting ROIs...')
  counts<-counts[which(counts$exon%in%rois$name),]  # reduce to target regions (from readCount)
  exons <- GRanges(seqnames = rois$chromosome, IRanges(start=rois$start,end=rois$end), names = rois$name)
  covered <- with(read.table(panel[1],header=FALSE),
                  GRanges(seqnames = V1, IRanges(start=V2+1, end=V3, names=V4),names=V4))  # BED FILE 0-based
  coveredexons<-subsetByOverlaps(exons,covered)
  ce<-data.frame(seqnames=seqnames(coveredexons), starts=start(coveredexons)-1, ends=end(coveredexons), name=mcols(coveredexons)$names)

  #
  # run exomedepth
  #
  # select normal samples (if 3+ specified)
  message('Selecting Normals (if any starting with NNN)...')
  refsamplenames<-as.vector(sapply(bam,basename))
  normals<-which(unlist(lapply(refsamplenames,function(x) substr(x,1,3)))=='NNN')
  if (length(normals)>2) refsamplenames<-refsamplenames[normals]

  #
  # Define low coverage exons
  #
  limit.coverage<-100
  exonnames<-coveredexons@elementMetadata@listData$names
  coverage.df<-data.frame(
    exon=counts$exon,
    gc=counts$GC,
    coverage.min=apply(counts[,refsamplenames],1,min),
    coverage.median=apply(counts[,refsamplenames],1,median),
    coverage.max=apply(counts[,refsamplenames],1,max))
  coverage.table<-coverage.df[which(coverage.df$coverage.median<limit.coverage & coverage.df$exon%in%exonnames),]

  #
  # get dispersion and correlation metrics
  #
  #refsamplestats<-lapply(refsamplenames, function(testsample) selectReferenceSet(counts, testsample, refsamplenames, TRUE)$summary.stats)
  #names(refsamplestats)<-refsamplenames

  # estimate reference set
  refsets[[samplename]]<-selectReferenceSet(counts, testsample, refsamplenames, TRUE)
  #
  # run exome depth
  #
  # prepare reference (sum reference choice)
  reference.selected<-apply(X = as.matrix(counts[, refsets[[samplename]]$reference.choice, drop = FALSE]),
                            MAR = 1, FUN = sum)
  reference.selected<-as.vector(reference.selected)
  # create ED object
  message('*** Creating ExomeDepth object...')
  suppressWarnings(ED <- new('ExomeDepth',
                             test = counts[,testsample],
                             reference = reference.selected,
                             formula = 'cbind(test, reference) ~ 1'))
  message('==========')
  # call CNV
  message('Calling CNVs...')
  result<-CallCNVs(x = ED,
                   transition.probability = threshold,
                   chromosome = counts$chromosome,
                   start = counts$start, 
                   end =counts$end,
                   name = counts$exon)

  # annotate results
  print(paste('raw CNV count:',length(result@CNV.calls)))

  # 
  message('Annotating CNVs...')
  if (length(result@CNV.calls)>0) {
    # add exon numbers (from subset)
    result.annotated<-AnnotateExtra(x = result, reference.annotation = coveredexons,
      min.overlap = 0.0001, column.name = 'exons.hg19')
    result.annotated@annotations$name<-as.factor(sapply(strsplit(as.character(result.annotated@annotations$name),'_'),"[[",1))
    results[[samplename]]<-result.annotated
    # write BED file of CNVs (is 1-based but irrelevant as only for visualisation)
    write.table(results[[samplename]]@CNV.calls[which(results[[samplename]]@CNV.calls[,9]>0),c(7,5,6,3,9)], file=sub("[.][^.]*$", ".bed", args[2], perl=TRUE), sep='\t', col.names=FALSE, quote=FALSE, row.names=FALSE)
  } else {
    # no results, empty bed file
    results[[samplename]]<-FALSE  # no CNVs called
    write.table(NULL, file=sub("[.][^.]*$", ".bed", args[2], perl=TRUE), sep='\t', quote=FALSE)
  }
} else {
  message('Skipping ExomeDepth (not enough reference samples)...')
  results[[samplename]]<-NA
  #
  # write BED file of CNVs (is 1-based but irrelevant as only for visualisation)
  #
  write.table(NULL, file=sub("[.][^.]*$", ".bed", args[2], perl=TRUE), sep='\t', quote=FALSE)
}
#
# knit report (using refsets, results)
#
knitrScript<-paste(scriptDirectory, 'exomeDepth.Rnw', sep='/')
print(knitrScript)
if (sub(".*[.]","",args[2],perl=TRUE)=="pdf") {
  knit2pdf(knitrScript, output=sub("[.][^.]*$", ".tex", args[2], perl=TRUE))
}
save.image(sub("[.][^.]*$", ".RData", args[2], perl=TRUE))
