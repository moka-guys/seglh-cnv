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
# load/combine exon/ROI data, unique and set to 1-based
#
rois<-NULL
for (tf in targetfiles) rois<-rbind(rois,read.table(tf,header=FALSE)[,1:4])
rois<-unique(rois)
rois[,2]<-rois[,2]+1
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
# save read count table as Rdata
# 
save(list=c("counts","rois","bam","bamstats"), file = args[1])

