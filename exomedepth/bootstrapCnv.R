# R script for exome depth data (takes a readcount data object created by readCount.R)
#
# LOAD LIBRARIES
#
cat(paste("Loading Packages\n"))
require(warn.conflicts=FALSE,quietly=FALSE,package="GenomicRanges")
require(warn.conflicts=FALSE,quietly=FALSE,package="ExomeDepth")
require(warn.conflicts=FALSE,quietly=FALSE,package="Rsamtools")
require(warn.conflicts=FALSE,quietly=TRUE,package="stringr")
require(warn.conflicts=FALSE,quietly=TRUE,package="dplyr")
require(warn.conflicts=FALSE,quietly=TRUE,package="moments")

#
# configs
#
sex <- regex("_[MF]_")
"%nin%" <- function(x, y) !("%in%"(x,y))
# shannon entropy to determine common CNVs
entropy <- function(target) {
  freq <- table(target) / length(target)
  vec <- as.data.frame(freq)[, 2]
  vec <- vec[vec > 0]
  -sum(vec * log2(vec))
}
#
# READ ARGS
#
args <- commandArgs(trailingOnly = TRUE)
setwd(dirname(args[1]))  # set working directory to target directory
load(args[1])
args <- commandArgs(trailingOnly = TRUE)
tries <- ifelse(!is.na(as.numeric(args[2])), as.numeric(args[2]), 100)
samplesize <- ifelse(!is.na(as.numeric(args[3])), as.numeric(args[3]), 30)
from <- ifelse(!is.na(as.numeric(args[4])), as.numeric(args[4]), 1)
to <- ifelse(!is.na(as.numeric(args[5])), as.numeric(args[5]),
  length(testsamplenames))

# load BED files to mark overlaps with known CNVs
beds <- args[which(endsWith(args, ".bed"))]
annotations <- NULL
if (length(beds) > 0) {
  message("Getting Segment annotations...")
  ranges <- list()
  for (bed in beds) {
    ranges[[bed]] <- with(
      read.table(bed, header = FALSE, fill = TRUE),
      GRanges(seqnames = V1, IRanges(start = V2 + 1, end = V3)))  # BED 0-based
  }
  annotations <- unlist(as(ranges, "GRangesList"))
}

#
# perform pre-filtering if scores files supplied
#
score_files <- args[which(endsWith(args, ".scores"))]
remove <- NULL
if (length(score_files) > 0) {
  proceed <- length(args[which(args == "RUN")]) != 0
  message("Analysing previous runs...")
  scores <- read.table(score_files[1], header = FALSE)
  if (length(score_files) > 1) {
    for (s in 2:length(score_files)) {
      scores <- rbind(scores, read.table(score_files[s]))
    }
  }
  colnames(scores) <- c("sample", "exon", "type", "prop", "common")
  # exclude common CNVs
  print(nrow(scores))
  scores <- scores[which(is.na(scores$common) | !scores$common), ]
  print(scores)
  # remove high-entropy, high-prevalence exons (likely common CNVs)
  scores %>% 
    group_by(exon) %>%
    summarise(entropy = entropy(prop), skewness = skewness(prop), n = n(),
      prevalence = n() / length(testsamplenames)) -> exon.entropy
  exon.entropy %>% print(n = Inf)
  exon.entropy %>% filter(entropy > 2 & prevalence > 0.4) %>% pull(exon) -> common
  print('Common CNVs in previous runs')
  print(common)
  scores <- scores[which(scores$exon %nin% common), ]
  # get samples with CNVs in more than half of bootstrap (likely true CNVs)
  scores <- scores[which(scores$prop > 0.5), ]
  exclusions <- unique(scores$sample)
  if (length(exclusions) == 0) {
    stop("Good set of normals, no further resampling required.")
  } else {
    print(scores)
    print(paste0("Excluding ", length(exclusions), "/", length(testsamplenames), ":"))
    print(exclusions)
    if (!proceed) stop("To proceed, run same command with extra 'RUN' argument")
    # remove samples and proceed as normally
    testsamplenames <- testsamplenames[which(testsamplenames %nin% exclusions)]
    refsamplenames$batch <- testsamplenames
  }
}

#
# read range of samples to bootstrap
#
from <- ifelse(!is.na(as.numeric(args[4])), as.numeric(args[4]), 1)
to <- ifelse(!is.na(as.numeric(args[5])), as.numeric(args[5]),
  length(testsamplenames))

#
# print what will be printed
#

message(paste("Read Count Data:", args[1]))

#
# Check if run in batch mode
#
if (!all(refsamplenames$batch == testsamplenames)) {
  stop(paste("You must provide data for batch normalisation."))
}

#
# run if enough samples
#
cnvcount <- vector()
refcount <- vector()
passed <- vector()
i <- 0
bed.data <- data.frame()
ref.data <- data.frame(
  test = c(),
  ref = c()
)

exons <- GRanges(seqnames = rois$chromosome,
  IRanges(start = rois$start, end = rois$end),
  names = rois$name)

timestamp <- format(Sys.time(), "%y%m%d%H%M%S")
cnv.file <- sub("[.][^.]*$", paste0(".", timestamp, ".scores"),
  args[1], perl = TRUE)
data.file <- sub("[.][^.]*$", paste0(".", timestamp, ".RData"),
  args[1], perl = TRUE)
for (testsample in testsamplenames[from:to]) {
  print(paste(i + from, to, sep = "/"))
  i <- i + 1
  # results
  samplecnvs <- data.frame()
  results <- data.frame(name = character(), exon = character(),
    type = character(), prop = numeric(), common = logical())
  #
  # get refsamples
  #
  refsamples <- refsamplenames$batch[which(refsamplenames$batch != testsample)]
  print(refsamplenames$batch)
  #
  # pick reference sample (sex-matched)
  #
  hasXY <- any(rois[, 1] == "X") || any(rois[, 1] == "Y")
  tssx <- str_extract(testsample, sex)
  rcsx <- str_extract(refsamples, sex)
  hasSex <- !(is.na(tssx) || any(is.na(rcsx)))
  if (hasXY && hasSex) {
    refsamples <- refsamples[which(rcsx == tssx)]
    print(paste("==> SEX_MATCHED SAMPLING", samplesize, "of", length(refsamples)))
  } else {
    print(paste("==> SAMPLING", samplesize, "of", length(refsamples)))
  }
  #
  # run boostrap
  #
  for (t in 1:tries) {
    print(paste0("** ", t, " of ", tries, " (", i, "/", to - from + 1, ")"))
    #
    # sample with replacement
    #
    refsamples.sample <- sample(refsamples, samplesize, replace = TRUE)
    refcounts <- as.matrix(counts[, refsamples.sample])
    refset <- suppressWarnings(select.reference.set(
      test.counts = counts[, testsample],
      reference.counts = refcounts,
      bin.length = (counts$end - counts$start)
    ))
    #
    # prepare reference (sum reference choice)
    #
    reference.selected <- apply(
      X = as.matrix(refcounts[, refset$reference.choice, drop = FALSE]),
      MAR = 1, FUN = sum)
    reference.selected <- as.vector(reference.selected)
    #
    # run exome depth
    #
    message("*** Creating ExomeDepth object...")
    suppressWarnings(ED <- new("ExomeDepth",
                              test = counts[,testsample],
                              reference = reference.selected,
                              formula = "cbind(test, reference) ~ 1"))
    #
    # call CNV
    #
    message("*** Calling CNVs...")
    result <- CallCNVs(x = ED,
                    transition.probability = 10^-4,
                    chromosome = counts$chromosome,
                    start = counts$start,
                    end = counts$end,
                    name = counts$exon)
    if (length(result@CNV.calls) > 0) {
      print(result@CNV.calls)
      cnvs <- result@CNV.calls[result@CNV.calls[,9]>0,]
      if (nrow(cnvs)>0) {
        samplecnvs <- rbind(samplecnvs, cnvs[, c(1, 2, 9, 10, 11, 12)])
      }
    }
  }

  print("************")
  print(samplecnvs)
  print("************")

  if (nrow(samplecnvs)>0) {
    dups <- table(unlist(apply(samplecnvs[which(samplecnvs$reads.ratio > 1),],
      1, function(x) x[1]:x[2])))/tries
    dels <- table(unlist(apply(samplecnvs[which(samplecnvs$reads.ratio < 1),],
      1, function(x) x[1]:x[2])))/tries
    for (roi in unique(c(names(dups), names(dels)))) {
      exon.name <- exons$names[as.numeric(roi)]
      cnvs.ovp <- ifelse(is.null(annotations), NA, 
        overlapsAny(exons[as.numeric(roi)], annotations))
      print(paste("** aggregating", roi))
      # aggregate data
      if (!is.na(dups[roi])) {
        dups.result <- data.frame(name = testsample, exon = exon.name,
          type = "duplication", prop = as.vector(dups[roi]), common = cnvs.ovp)
        results <- rbind(results, dups.result)
      }
      if (!is.na(dels[roi])) {
        dels.result <- data.frame(name = testsample, exon = exon.name,
          type = "deletion", prop = as.vector(dels[roi]), common = cnvs.ovp)
        results <- rbind(results,dels.result)
      }
    }
    write.table(results, file=cnv.file, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE, append=TRUE)
    save.image(data.file)
  }
}