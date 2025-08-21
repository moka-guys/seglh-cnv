# R script for exome depth data (takes a readcount data object created by readCount.R)

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
require(warn.conflicts=FALSE,quietly=FALSE,package="dplyr")
require(warn.conflicts=FALSE,quietly=FALSE,package="tidyr")
source('ed2vcf.R')

#
# READ ARGS
#
cmd <- commandArgs(trailingOnly = FALSE)
runningScript <- unlist(strsplit(cmd[which(substr(cmd,1, 7) == "--file=")],"="))[2]
scriptDirectory <- normalizePath(dirname(runningScript))
args <- commandArgs(trailingOnly = TRUE)
setwd(dirname(args[2]))  # set working dir to target dir (ABSOLUTELY NECESSARY)

#
# print what will be printed
#
panel <- unlist(strsplit(args[3],":"))  # covered regons name
testsample <- basename(unlist(strsplit(args[5],":"))[1])  # RGID
samplename <- as.character(unlist(strsplit(args[5],":"))[2]) # RGSM
threshold <- as.numeric(unlist(strsplit(args[5],":"))[3]) # Thresh in ngs_config
if (is.na(threshold)) {
    threshold <- 10^-3
}

# version report header
argversion <- args[1]
versionfile <- paste(scriptDirectory, "VERSION", sep = "/")
pipeversion <- ifelse(file.exists(versionfile),
                    paste(argversion, readChar(versionfile, 7), sep = "-"),
                    argversion)
message(paste("Running", pipeversion))

# loads counts,samplenames,rois (extended exons.hg19)
load(args[4])
extras<-args[6]

# get bed file for filtering
df_bed <- read.table(panel[1], header = FALSE, sep = "\t", col.names = c("chromosome", "start", "end", "name"))
df_bed$chromosome <- as.character(df_bed$chromosome)
rois$chromosome <- as.character(rois$chromosome)
# merged cnv bed with rois that is same as readcount bed file in order to get the exon number info 
# which is present only in readcount bed but not in CNV bed
df_merged <- left_join(df_bed, rois, by = c("chromosome", "start", "end"))
# check if any regions present in CNV bed are absent in readcount bed
# if any regions in cnv bed are absent in readcount bed, replace NA with Unknown_gene
if(any(is.na(df_merged$name.y))) {
  print("WARNING: some regions in CNV bed are absent in readcount bed")
  print(df_merged[is.na(df_merged$name.y), ])
  df_merged$name.y[is.na(df_merged$name.y)] <- "Unknown_gene"
  }

# rename cols
df_merged <- df_merged %>%
                rename(
                gene = name.x,
                gene_exonnum = name.y
              )
# remove duplicated rows
df_bed_nodup <- df_merged[!duplicated(df_merged), ]

#
# report run configuration
#
message(paste('        Version:',args[1]))
message(paste('         Output:',args[2]))
message(paste('            ROI:',panel[1]))
message(paste('     Panel Name:',panel[2]))
message(paste('    Read Counts:',args[4]))
message(paste('   Original BAM:',testsample))
message(paste('     SampleName:',samplename))
message(paste('  Normalisation:',names(refsamplenames)))
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
if (length(names(refsamplenames))>0) {
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
  selected.genes<-unique(sapply(strsplit(df_bed_nodup$gene_exonnum,'_'),function(x) x[1]))

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
  annotations<-NA  # CNV annotations
  ## load from file (including annotations)
  if (!is.na(extras)) {
    extras<-ifelse(startsWith(extras,'/'), extras, paste(scriptDirectory, extras, sep='/'))
    load(extras)
  }

  #
  # Define low coverage exons
  #
  message('Identifying low coverage exons...')
  exonnames<-coveredexons@elementMetadata@listData$names
  coverage.df<-data.frame()
  for (rs in names(refsamplenames)) {
    coverage.df<-rbind(coverage.df,data.frame(
      refset=rs,
      exon=counts$exon,
      gc=counts$GC,
      coverage.min=apply(counts[,refsamplenames[[rs]]],1,min),
      coverage.median=apply(counts[,refsamplenames[[rs]]],1,median),
      coverage.max=apply(counts[,refsamplenames[[rs]]],1,max)
    ))
  }
  coverage.table<-coverage.df[which(coverage.df$coverage.median<limits$coverage & coverage.df$exon%in%exonnames),]

  #
  # check if reference sample choice has mismatched sex
  #
  message('Checking sex match...')
  sexMismatch<-list()
  for (rs in names(refsamplenames)) {
    sexMismatch[[rs]]<-rep(NA,length(refsets[[testsample]][[rs]]$reference.choice))
    if (any(rois[,1]=="X") || any(rois[,1]=="Y")) {
      sex<-regex('_[FM]_')
      # get sex of testsample (ignore undefined)
      tssx<-str_extract(testsample,sex)
      rcsx<-str_extract(refsets[[testsample]][[rs]]$reference.choice,sex)
      # if sex information available
      if (!is.na(tssx) && !all(is.na(rcsx))) {
        sexMismatch[[rs]]<-!(rcsx==tssx)
      }
    }
  }

  #
  # prepare reference (sum reference choice)
  #
  message('Preparing reference set...')
  reference.selected<-list()
  ref.correlation<-list()
  for (rs in names(refsamplenames)) {
    reference.selected[[rs]]<-as.vector(apply(
      X = as.matrix(counts[, refsets[[testsample]][[rs]]$reference.choice, drop = FALSE]),
      MAR = 1,
      FUN = sum
    ))
    ref.correlation[[rs]]<-cor(cbind(rpkm[,testsample],calcRPKM(reference.selected[[rs]],(counts$end-counts$start+1))))[1,2]
  }

  #
  # Build QC table
  #
  message('Building QC summary...')
  # decides if a particular metric is within the reference range
  # v being the measured value, t[1] being the warning limit, t[2] being the failure limit
  decide<-function(v,t) {
    result<-list(
      value=round(v,3),
      threshold='no threshold set',
      warning=NA,
      fail=NA,
      status="PASS"
    )
    if (!is.null(t)) {
      # comparison function (infers if lower or upper limit)
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
  save.image()
  message('•••SAVED•••')

  # create a qc table for each refset
  qc <- list()
  for (rs in unique(stats$refset)) {
    stats.refset<-stats[which(stats$refset == rs),
      which(colnames(stats) != "refset")]

    # calc CV Z-score for each batch
    batches<-as.factor(unlist(lapply(strsplit(rownames(stats.refset),'_'),'[[',1)))
    batches.mean<-by(stats.refset$coeff.var, batches, mean)
    batches.sd<-by(stats.refset$coeff.var, batches, sd)
    
    # get/calc testsample data
    testsample.cv<-stats.refset[which(stats.refset$sample==testsample),'coeff.var']
    testsample.batch<-batches[which(stats.refset$sample==testsample)]
    batch.mean<-batches.mean[[testsample.batch]]
    batch.sd<-batches.sd[[testsample.batch]]
    testsample.cvz<-(testsample.cv - batch.mean)/batch.sd

    qc[[rs]]<-rbind(
      decide(stats.refset[testsample,"batch.mediancor"],limits$medcor),
      decide(stats.refset[testsample,"batch.maxcor"],limits$maxcor),
      decide(stats.refset[testsample,"coeff.var"],limits$coeffvar),
      decide(testsample.cvz,limits$cvzscore),
      decide(ref.correlation[[rs]], limits$refcor),
      decide(stats.refset[testsample,"refsamples"],limits$refcount),
      decide(stats.refset[testsample,"min.refs"],limits$minrefs),
      decide(stats.refset[testsample,"expected.BF"],limits$expectedbf)
    )
    rownames(qc[[rs]])=c(
              "Median correlation in batch",
              "Maximum correlation in batch",
              "Coefficient of variation (CV)",
              "CV Z-score within batch",
              "Correlation with reference",
              "Size of reference set",
              "Forced minimum reference set size",
              "Expected BF with reference"
    )
    qc[[rs]]<-as.data.frame(qc[[rs]])
    print(qc[[rs]])
  }

  #
  # run exome depth
  #
  cnvs_all<-list()
  for (rs in names(reference.selected)) {
    message('*** Creating ExomeDepth object...')
    suppressWarnings(ED <- new('ExomeDepth',
                               test = counts[,testsample],
                               reference = reference.selected[[rs]],
                               formula = 'cbind(test, reference) ~ 1'))

    # call CNV
    message('*** Calling CNVs...')
    cnvs<-CallCNVs(x = ED,
                    transition.probability = threshold,
                    chromosome = counts$chromosome,
                    start = counts$start, 
                    end =counts$end,
                    name = counts$exon)
    print(paste('Raw CNV count:',length(cnvs@CNV.calls)))

    # annotate results
    message("*** Annotating CNVs...")
    exon_overlap_frac <- 0.000000001  # 1bp/1Gb, basically any overlap
    extra_overlap_frac <- 0.1  # 1bp/10bp (only valid for if annotation is large known SegDups)
    if (nrow(cnvs@CNV.calls) > 0) {
      message(paste("CNVs to be annotated:", length(cnvs@CNV.calls)))
      # add exon numbers (from subset)
      cnvs.annotated <- AnnotateExtra(x = cnvs,
                                    reference.annotation = coveredexons,
                                    min.overlap = exon_overlap_frac,
                                    column.name = "exons.hg19")

      # add extra annotation
      if (all(is.na(annotations))) {
        message("No CNV annotaions available")
      } else {
        message("Adding extra annotations...")
        cnvs.annotated <- AnnotateExtra(x = cnvs.annotated,
                                        reference.annotation = annotations,
                                        min.overlap = extra_overlap_frac,
                                        column.name = "annotation")
      }
      cnvs.annotated@annotations$name <- as.factor(
        sapply(strsplit(
          as.character(cnvs.annotated@annotations$name), '_'), "[[", 1))
      
      if (nrow(cnvs.annotated@CNV.calls) > 0) {
        # get the non NA row which are actual CNV call 
        cnvcall <- cnvs.annotated@CNV.calls[!is.na(cnvs.annotated@CNV.calls$exons.hg19), ]
        if (nrow(cnvcall) > 0) {
          updated_cnv <- data.frame()
          for (i in 1:nrow(cnvcall)) {
            # subset each row
            df_subset <- cnvcall[i:i,]
            # expand the row for individual exon number
            cnvcall_expanded <- df_subset %>%
                              separate_rows(exons.hg19, sep = ",")
            # join with bed file to filter out regions that are not present in cnv bed
            cnv_bed <- inner_join(cnvcall_expanded, df_bed_nodup, by = c("exons.hg19" = "gene_exonnum"))
            
            # replace the start and end with correct interval
            cnv_bed$start.x <-min(cnv_bed$start.y)+1 # to make 1 based interval
            cnv_bed$end.x <-max(cnv_bed$end.y)
            # drop unnecessary col
            cnv_bed <- subset(cnv_bed, select = -c(id, chromosome.y, start.y, end.y, gene))
            # replace id col with correct interval
            cnv_bed$id <- paste0("chr",cnv_bed$chromosome.x,":", cnv_bed$start.x,"-", cnv_bed$end.x)
            # rename the col
            cnv_bed <- cnv_bed %>%
                    rename(
                    start = start.x,
                    end = end.x,
                    chromosome = chromosome.x
                  )
            # collapse back into single row
            cnv_summary <- cnv_bed %>%
                        group_by(across(-exons.hg19)) %>%
                        summarise(exons.hg19= paste(exons.hg19, collapse = ","), .groups = "drop")
            # append into updated_csv
            updated_cnv <- rbind(updated_cnv, cnv_summary)
          }
          # replace old cnv call regions with updated cnv regions
          cnvs.annotated@CNV.calls <- cnvs.annotated@CNV.calls %>%
                                 filter(!exons.hg19 %in% updated_cnv$exons.hg19) %>%
                                 bind_rows(updated_cnv)
        }
      }

      results[[rs]] <- cnvs.annotated
    } else {
      message(paste('No CNVs called with refset',toupper(rs)))
      results[[rs]] <- NA  # no CNVs called
    }

    # write BED file of CNVs
    message(paste("Writing BED file for refset", toupper(rs)))
    getBed <- function(x) {
      if (class(x) == "ExomeDepth") {
        beddata <- x@CNV.calls[which(x@CNV.calls[, 9] > 0), c(7, 5, 6, 3, 9)]
        beddata[, 2] <- beddata[, 2] - 1
        return(beddata)
      }
      NULL
    }
    bed.filename <- sub("[.][^.]*$", paste0(".", rs, ".bed"), args[2],
      perl = TRUE)
    write.table(
      getBed(results[[rs]]),
      file = bed.filename, sep = "\t",
      col.names = FALSE, quote = FALSE, row.names = FALSE)

    # write VCF file
    message(paste("Writing VCF file for refset", toupper(rs)))
    if (file.exists(referenceFasta)) {
      vcf.filename <- sub("[.][^.]*$", paste0(".", rs, ".vcf"), args[2],
        perl=TRUE)
      ed2vcf(results[[rs]], vcf.filename, referenceFasta, rois, samplename)
    }
  }
} else {
  message("Skipping ExomeDepth. No suitable reference set available")
}

#
# knit report (using refsets, results)
#
knitrScript <- paste(scriptDirectory, "exomeDepth.Rnw", sep="/")
if (sub(".*[.]","",args[2],perl=TRUE)=="pdf") {
  knit2pdf(knitrScript, output=sub("[.][^.]*$", ".tex", args[2], perl=TRUE))
}
save.image(sub("[.][^.]*$", ".RData", args[2], perl = TRUE))
