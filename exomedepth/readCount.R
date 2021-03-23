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
#
# READ ARGS
#
args<-commandArgs(trailingOnly = TRUE)

#
# print what will be printed
#
message(paste('         Output:',args[1]))
message(paste('      Reference:',args[2]))
# read ROI
targetfiles<-vector()
for (roi in unlist(strsplit(args[3],','))) {
    targetfiles<-append(targetfiles,roi)
    message(paste('            ROI:',roi))
}

# read tracks (BAM files)
bam<-args[4:length(args)]
for (bamfile in bam) {
    message(paste('          Track:',bamfile))
}

#
# load/combine exon/ROI data and make unique
#
rois<-NULL
for (tf in targetfiles) rois<-rbind(rois,read.table(tf,header=FALSE)[,1:4])
rois<-unique(rois)
colnames(rois)<-c('chromosome','start','end','name')

#
# generate read counts from BAM files for all exons
#
suppressWarnings(counts <- getBamCounts(
                            bed.frame = rois,
                            bam.files = bam,
                            min.mapq = 10,
                            include.chr = FALSE,  # chrom start with chr prefix
                            referenceFasta = args[2]))

counts<-as(counts, 'data.frame')

#
# Calculate RPKM and batch statistics
#
calcRPKM<-function(c,l) c/(l*sum(c)/10^6)
counts.len<-counts$end-counts$start+1
rpkm<-apply(counts[,c(6:ncol(counts))],2,function(x) calcRPKM(x,counts.len))
batch.cv<-apply(rpkm,2,function(r) sd(r)/mean(r)*100)
batch.cor<-cor(rpkm)
diag(batch.cor)<-NA

#
# Generate sample correlation/dispersion data
#
pickreference<-function(testsample) {
  refsamples<-bamnames[which(bamnames!=testsample)]
  select.reference.set(
    test.counts=counts[,testsample],
    reference.counts=as.matrix(counts[,refsamples]),
    bin.length=(counts$end - counts$start)
  )$summary.stats
}
bamnames<-as.vector(sapply(bam,basename))
bamstats<-lapply(bamnames, pickreference)
names(bamstats)<-bamnames

#
# statistics output
#
stats<-data.frame()
for (testsample in colnames(rpkm)) {
    d<-cbind(
             sample=testsample,
             refsamples=which(bamstats[[testsample]]$selected),
             bamstats[[testsample]][which(bamstats[[testsample]]$selected),
                                    c('correlations','expected.BF','phi','RatioSd','mean.p','median.depth')],
             batch.maxcor=max(batch.cor[testsample,],na.rm=TRUE),
             batch.mediancor=median(batch.cor[testsample,],na.rm=TRUE),
             coeff.var=batch.cv[testsample]
            )
    rownames(d)<-testsample
    stats<-rbind(stats,d)
}
write.table(stats, file=sub("[.][^.]*$", ".csv", args[1], perl=TRUE), sep='\t', quote=FALSE, row.names=FALSE)

#
# save read count table as Rdata
# 
save(list=c("counts","rois","bam","bamstats","batch.cv","batch.cor","rpkm","calcRPKM","stats"), 
     file = args[1])

