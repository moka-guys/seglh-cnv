# R script for exome depth data
#
# INSTALL PACKAGES
#
options(menu.graphics=FALSE)
options(repos=structure(c(CRAN="http://cran.ma.imperial.ac.uk/")))
if (!("ExomeDepth" %in% installed.packages())) {
  install.packages("ExomeDepth")
}
#
# LOAD LIBRARIES
#
cat(paste('Loading Packages\n'))
require(warn.conflicts=FALSE,quietly=TRUE,package="plyr")
require(warn.conflicts=FALSE,quietly=TRUE,package="ExomeDepth")
require(warn.conflicts=FALSE,quietly=TRUE,package="stringr")

#
# configs
#
normalprefix<-'NORMAL'
cmp.cols<-3
sex<-regex('_[MF]_')

#
# READ ARGS
#
args<-commandArgs(trailingOnly = TRUE)

#
# print what will be printed
#
message(paste('         Output:',args[1]))
message(paste('      Reference:',args[2]))
referenceFasta<-args[2]
# read ROI
targetfiles<-vector()
for (roi in unlist(strsplit(args[3],','))) {
    targetfiles<-append(targetfiles,roi)
    message(paste('            ROI:',roi))
}

#
# load/combine exon/ROI data and make unique
#
rois<-NULL
for (tf in targetfiles) rois<-rbind(rois,read.table(tf,header=FALSE)[,1:4])
rois<-unique(rois)
colnames(rois)<-c('chromosome','start','end','name')


# read tracks (BAM files, precomputed normals)
bam<-vector()
pcn<-vector()
for (bamfile in args[4:length(args)]) {
    if (endsWith(bamfile,'.bam')) {
        bam<-append(bam,bamfile)
        message(paste('     Alignments:',bamfile))
    } else if (endsWith(bamfile,'.RData')) {
        pcn<-append(pcn,bamfile)
        message(paste('         Counts:',bamfile))
    } else {
        message(paste('        Skipped:',bamfile))
    }
}

#
# omit normals that are in supplied counts file
#
names.imported<-vector()
for (c in pcn) {
    attach(c)
    names.imported<-append(names.imported,
                           colnames(counts)[which(startsWith(colnames(counts),normalprefix))])
    detach()
}
# remove imported normals from bam list
bams.to.count<-which(!as.vector(sapply(bam,basename))%in%names.imported)
bam<-bam[bams.to.count]

#
# generate read counts from BAM files for all exons
#
suppressWarnings(counts.computed <- getBamCounts(
                            bed.frame = rois,
                            bam.files = bam,
                            min.mapq = 10,
                            include.chr = FALSE,  # chrom start with chr prefix
                            referenceFasta = referenceFasta))
counts.computed<-as(counts.computed, 'data.frame')


#
# add precomputed counts (normals only)
# 
for (c in pcn) {
    # import precomputed normals
    attach(c)
    targets.imported<-counts[,c(1:cmp.cols)]
    counts.imported<-counts[,which(startsWith(colnames(counts),normalprefix))]
    detach()
    # if compatible normals, amend counts table
    targets.same<-all.equal(counts.computed[,c(1:cmp.cols)], targets.imported)
    if (targets.same && ncol(counts.imported)>0) {
        counts.computed<-cbind(counts.computed,counts.imported)
    } else {
        message(paste('ERROR: Supplied count data not compatible',c))
    }
}
counts<-counts.computed

#
# select normal samples (if 3+ specified)
#
refsamplenames<-colnames(counts)[which(endsWith(colnames(counts),'.bam'))]
testsamplenames<-refsamplenames
# NORMAL selection
normals<-which(substr(refsamplenames,1,nchar(normalprefix)) == normalprefix)
tests<-which(substr(testsamplenames,1,nchar(normalprefix)) != normalprefix)
if (length(normals)>0) {
    message(paste("-> Panel of Normals [",length(normals),"]..."))
    refsamplenames<-refsamplenames[normals]
    testsamplenames<-testsamplenames[tests]
    # ensure there are 0 (pon), 2 or more (batch) samples
    if (length(testsamplenames)==1) {
        message('NOTE: Single test sample. Batch correlation will be 1.')
        testsamplenames<-rep(testsamplenames,2)
    }
} else {
    message('-> Batch normalisation')
}

#
# Process test samples if any
#
if (length(testsamplenames)==0) {
    message('Creating Panel of Normals (PoN)...')
    # save normal counts
    save(list=c(
                "refsamplenames",  # supplied normals
                "counts"),         # normal counts
         file = args[1])
} else {
    message('Building reference sets...')
    # Calculate RPKM
    calcRPKM<-function(c,l) c/(l*sum(c)/10^6)
    counts.len<-counts$end-counts$start+1

    rpkm<-apply(counts[,testsamplenames],2,function(x) calcRPKM(x,counts.len))

    # Calculate batch statistics (all samples excluding normals if any)
    batch.cv<-apply(rpkm,2,function(r) sd(r)/mean(r)*100)
    batch.cor<-cor(rpkm)
    diag(batch.cor)<-NA

    # Pick reference sample set
    selectReferenceSet<-function(testsample) {
      message(paste('Picking reference for',testsample))
      # exclude test sample
      refsamples<-refsamplenames[which(refsamplenames!=testsample)]
      # if PoN && XY && Sex -> sex match
      hasPon<-length(normals)>0
      hasXY<-any(rois[,1]=="X") || any(rois[,1]=="Y")
      tssx<-str_extract(testsample,sex)
      rcsx<-str_extract(refsamples,sex)
      hasSex<-!(is.na(tssx) || any(is.na(rcsx)))
      if (hasPon && hasXY && hasSex) {
        message(paste('Enforcing sex match of', testsample, 'to reference'))
        refsamples<-refsamples[which(rcsx==tssx)]
      }
      # pick reference  
      suppressWarnings(select.reference.set(
        test.counts=counts[,testsample],
        reference.counts=as.matrix(counts[,refsamples]),
        bin.length=(counts$end - counts$start)
      ))
    }
    refsets<-lapply(testsamplenames, selectReferenceSet)
    names(refsets)<-testsamplenames

    # statistics output
    stats<-data.frame()
    for (testsample in testsamplenames) {
        message(paste('Collecting statistics for',testsample))
        d<-cbind(
                 sample=testsample,
                 refsamples=which(refsets[[testsample]]$summary.stats$selected),
                 refsets[[testsample]]$summary.stats[which(refsets[[testsample]]$summary.stats$selected),
                                        c('correlations','expected.BF','phi','RatioSd','mean.p','median.depth')],
                 batch.maxcor=max(batch.cor[testsample,],na.rm=TRUE),
                 batch.mediancor=median(batch.cor[testsample,],na.rm=TRUE),
                 coeff.var=batch.cv[testsample]
                )
        rownames(d)<-testsample
        stats<-rbind(stats,d)
    }
    print(stats)
    write.table(stats, file=sub("[.][^.]*$", ".csv", args[1], perl=TRUE), sep='\t', quote=FALSE, row.names=FALSE)

    # save read count table as Rdata
    save(list=c(
                "referenceFasta",  # reference fasta
                "refsamplenames",  # selected reference samples (normals if provided)
                "testsamplenames", # supplied test samples
                "counts",          # all read counts
                "rois",            # regions of interest (full)
                "batch.cv",        # Coefficient of Variation for each sample in batch
                "batch.cor",       # Correlation within batch
                "refsets",         # picked reference sets
                "rpkm",            # RPKM calculation
                "calcRPKM",        # function to calculate RPKM
                "stats"),          # model and betch statistics
         file = args[1])

}
message('DONE')
