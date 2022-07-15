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
require(warn.conflicts=FALSE,quietly=FALSE,package="xtable")
require(warn.conflicts=FALSE,quietly=FALSE,package="knitr")
require(warn.conflicts=FALSE,quietly=FALSE,package="kableExtra")
require(warn.conflicts=FALSE,quietly=FALSE,package="randomForest")
require(warn.conflicts=FALSE,quietly=FALSE,package="Rsamtools")
require(warn.conflicts=FALSE,quietly=FALSE,package="stringr")
source('ed2vcf.R')

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
    threshold<-10^-3
}

# version report header
argversion<-args[1]
versionfile<-paste(scriptDirectory,'VERSION', sep='/')
pipeversion<-ifelse(file.exists(versionfile),
                    paste(argversion,readChar(versionfile,7),sep="-"),
                    argversion)
message(paste('Running',pipeversion))

# loads counts,samplenames,rois (extended exons.hg19)
load(args[4])
normalisation.method<-ifelse(testsample%in%refsamplenames,'BATCH','PoN')
extras<-args[6]
message(paste('        Version:',args[1]))
message(paste('         Output:',args[2]))
message(paste('            ROI:',panel[1]))
message(paste('     Panel Name:',panel[2]))
message(paste('    Read Counts:',args[4]))
message(paste('   Original BAM:',testsample))
message(paste('     SampleName:',samplename))
message(paste('  Normalisation:',normalisation.method))
message(paste('    Ref samples:',length(refsamplenames[which(testsample!=refsamplenames)])))
message(paste('         Extras:',extras))

#
# Check if testsample in refsamples
#
if (!testsample%in%names(refsets)) {
  stop(paste("The requested sample not available in", paste(names(refsets),collapse=',')))
}

#
# run if enough samples
#
results<-list()
qc<-data.frame() ## to make sure qc table gets saved
if (length(refsamplenames)>=3) {
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
  selected.genes<-unique(sapply(strsplit(ce$name,'_'),function(x) x[1]))

  #
  # set defaults and load QC limits (override)
  #
  limits<-list(
    medcor=c(NA, 0.90),    # median correlation within batch
    maxcor=c(0.95, 0.90),  # max correlation within batch
    refcor=c(0.95, 0.90),  # reference set correlation
    refcount=c(3,1),       # refernce set size (selected reference samples)
    coeffvar=c(NA, NA),    # coefficient of variation
    cvzscore=c(1,4),       # Z-Score of the coefficient of variation 
    coverage=c(100),       # Minimum exon depth (read count)
    expectedbf=c(5.0, NA), # expected BF
    minrefs=c(2,Inf)       # minimum reference set size
  )
  predicted_qc<-NA
  annotations<-NA  # CNV annotations
  ## load from file
  if (!is.na(extras)) {
    extras<-ifelse(startsWith(extras,'/'), extras, paste(scriptDirectory, extras, sep='/'))
    load(extras)
  }


  #
  # Define low coverage exons
  #
  exonnames<-coveredexons@elementMetadata@listData$names
  coverage.df<-data.frame(
    exon=counts$exon,
    gc=counts$GC,
    coverage.min=apply(counts[,refsamplenames],1,min),
    coverage.median=apply(counts[,refsamplenames],1,median),
    coverage.max=apply(counts[,refsamplenames],1,max))
  coverage.table<-coverage.df[which(coverage.df$coverage.median<limits$coverage & coverage.df$exon%in%exonnames),]

  #
  # check if reference sample choice has mismatched sex
  #
  sexMismatch<-rep(NA,length(refsets[[testsample]]$reference.choice))
  if (any(rois[,1]=="X") || any(rois[,1]=="Y")) {
    sex<-regex('_[FM]_')
    # get sex of testsample (ignore undefined)
    tssx<-str_extract(testsample,sex)
    rcsx<-str_extract(refsets[[testsample]]$reference.choice,sex)
    # if sex information available
    if (!is.na(tssx) && !all(is.na(rcsx))) {
      sexMismatch<-!(rcsx==tssx)
    }
  }
  #
  # prepare reference (sum reference choice)
  #
  reference.selected<-apply(X = as.matrix(counts[, refsets[[testsample]]$reference.choice, drop = FALSE]),
                            MAR = 1,
                            FUN = sum)
  reference.selected<-as.vector(reference.selected)

  #
  # Build QC table and predict Quality outcome of classifier provided
  #
  if (!is.na(extras) && !is.null(rfc)) {
    predicted_qc<-predict(rfc,stats[testsample,2:ncol(stats)])
  }
  # build QC table
  ref.correlation<-cor(cbind(rpkm[,testsample],calcRPKM(reference.selected,(counts$end-counts$start+1))))[1,2]
  # decide on failures
  decide<-function(v,t) {
    result<-list(
      value=round(v,3),
      threshold='no threshold set',
      warning=NA,
      fail=NA,
      status="PASS"
    )
    if (!is.null(t)) {
      cmp<-ifelse(any(is.na(t)) || t[1]>t[2], function(m,n) m>=n, function(m,n) m<n)
      threshold<-ifelse(any(is.na(t)) || t[1]>t[2],"equal or greater than","less than" )
      status<-ifelse(!is.na(t[2]) && !cmp(v,t[2]),"FAIL",
                  ifelse(!is.na(t[1]) && !cmp(v,t[1]),"CAUTION","PASS"))
      result<-list(
        value=round(v,3),
        threshold=threshold,
        warning=t[1],
        fail=t[2],
        status=status
      )
    }
    result
  }

  # calc CV Z-score for each batch
  batches<-as.factor(unlist(lapply(strsplit(rownames(stats),'_'),'[[',1)))
  batches.mean<-by(stats$coeff.var, batches, mean)
  batches.sd<-by(stats$coeff.var, batches, sd)
  
  # get/calc testsample data
  testsample.cv<-stats[which(rownames(stats)==testsample),'coeff.var']
  testsample.batch<-batches[which(rownames(stats)==testsample)]
  batch.mean<-batches.mean[[testsample.batch]]
  batch.sd<-batches.sd[[testsample.batch]]
  testsample.cvz<-(testsample.cv - batch.mean)/batch.sd

  # get Z-score bands
  cvz.bands<-list(
    c(batch.mean-(4*batch.sd),batch.mean+(4*batch.sd)), 
    c(batch.mean-(3*batch.sd),batch.mean+(3*batch.sd)), 
    c(batch.mean-(2*batch.sd),batch.mean+(2*batch.sd)), 
    c(batch.mean-(1*batch.sd),batch.mean+(1*batch.sd))
  )

  qc<-rbind(
    decide(stats[testsample,"batch.mediancor"],limits$medcor),
    decide(stats[testsample,"batch.maxcor"],limits$maxcor),
    decide(stats[testsample,"coeff.var"],limits$coeffvar),
    decide(testsample.cvz,limits$cvzscore),
    decide(ref.correlation, limits$refcor),
    decide(stats[testsample,"refsamples"],limits$refcount),
    decide(stats[testsample,"min.refs"], limits$minrefs),
    decide(stats[testsample,"expected.BF"],limits$expectedbf)
  )
  rownames(qc)=c(
            "Median correlation in batch",
            "Maximum correlation in batch",
            "Coefficient of variation (CV)",
            "CV Z-score within batch",
            "Correlation with reference",
            "Size of reference set",
            "Forced minimum reference set size",
            "Expected BF with reference"
  )
  qc<-as.data.frame(qc)
  print(qc)

  #
  # run exome depth
  #
  message('*** Creating ExomeDepth object...')
  suppressWarnings(ED <- new('ExomeDepth',
                             test = counts[,testsample],
                             reference = reference.selected,
                             formula = 'cbind(test, reference) ~ 1'))
  # call CNV
  message('*** Calling CNVs...')
  result<-CallCNVs(x = ED,
                   transition.probability = threshold,
                   chromosome = counts$chromosome,
                   start = counts$start, 
                   end =counts$end,
                   name = counts$exon)

  # annotate results
  print(paste('Raw CNV count:',length(result@CNV.calls)))
  message('Annotating CNVs...')
  exon_overlap_frac<-0.000000001  # 1bp/1Gb, basically any overlap
  extra_overlap_frac<-0.1  # 1bp/10bp (only valid for if annotation is large known SegDups)
  if (length(result@CNV.calls)>0) {
    # add exon numbers (from subset)
    result.annotated<-AnnotateExtra(x = result, reference.annotation = coveredexons,
      min.overlap = exon_overlap_frac, column.name = 'exons.hg19')
    # add extra annotation
    if (all(is.na(annotations))) {
      message('No CNV annotaions available')
    } else {
      message('Adding extra annotations...')
      result.annotated<-AnnotateExtra(x = result.annotated, reference.annotation = annotations,
                                      min.overlap = extra_overlap_frac, column.name = 'annotation')
    }
    result.annotated@annotations$name<-as.factor(sapply(strsplit(as.character(result.annotated@annotations$name),'_'),"[[",1))
    results[[testsample]]<-result.annotated
    # write BED file of CNVs
    bed.data<-results[[testsample]]@CNV.calls[which(results[[testsample]]@CNV.calls[,9]>0),c(7,5,6,3,9)]
    bed.data[,2]<-bed.data[,2]-1
    write.table(bed.data, file=sub("[.][^.]*$", ".bed", args[2], perl=TRUE), sep='\t', col.names=FALSE, quote=FALSE, row.names=FALSE)
    # write VCF file
    if (file.exists(referenceFasta)) {
      ed2vcf(results[[testsample]], sub("[.][^.]*$",".vcf",args[2],perl=TRUE), referenceFasta, rois, samplename)
    }
  } else {
    # no results, empty BED and VCF files
    results[[testsample]]<-FALSE  # no CNVs called
    write.table(NULL, file=sub("[.][^.]*$", ".bed", args[2], perl=TRUE), sep='\t', quote=FALSE)
    if (file.exists(referenceFasta)) {
      ed2vcf(NULL, sub("[.][^.]*$",".vcf",args[2],perl=TRUE), referenceFasta, rois, samplename)
    }
  }
} else {
  message('Skipping ExomeDepth (not enough reference samples)...')
  results[[testsample]]<-NA
}

#
# knit report (using refsets, results)
#
knitrScript<-paste(scriptDirectory, 'exomeDepth.Rnw', sep='/')
if (sub(".*[.]","",args[2],perl=TRUE)=="pdf") {
  knit2pdf(knitrScript, output=sub("[.][^.]*$", ".tex", args[2], perl=TRUE))
}
save.image(sub("[.][^.]*$", ".RData", args[2], perl=TRUE))
