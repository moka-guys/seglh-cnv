require(Rsamtools)
require(GenomicRanges)

ed2vcf<-function(x, filename, fasta, roi, samplename) {
    # Modelled after: https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/CNVVCFFile_fDG_dtSW.htm
    ## Define Header
    header<-list(
      fileformat="VCFv4.2",
      fileDate = format(Sys.time(), "%Y-%m-%d"),
      source = "exomedepth",
      reference = fasta,
      ALT = data.frame(
        ID=c("CNV","DEL","DUP"),
        Description=c(
          "Copy number variant region",
          "Deletion relative to the reference",
          "Region of elevated copy number relative to the reference")
      ),
      FILTER = data.frame(
        ID=c("outOfScope","cnvLowBF","cnvCopyRatio"),
        Description=c(
          "CNV out of scope",
          "CNV with BF below 3",
          "CNV with copy ratio within +/- 0.25 of 1.0")
      ),
      INFO = data.frame(
        ID=c("SVLEN","SVTYPE","END","CIPOS","CIEND","GENES"),
        Number=c(".","1","1","2","2","."),
        Type=c('Integer','String','Integer','Integer','Integer','String'),
        Description=c(
          "Difference in length between REF and ALT alleles",
          "Type of structural variant",
          "End position of the variant described in this record",
          "Confidence interval around POS",
          "Confidence interval around END",
          "Affected genes (regions)")
      ),
      FORMAT=data.frame(
        ID=c("GT","SM","CN","BC","BF","EL"),
        Number=c("1","1","1","1","1","1"),
        Type=c("String","Float","Integer","Integer","Float","String"),
        Description=c(
          "Genotype",
          "Linear copy ratio of the segment mean",
          "Estimated copy number",
          "Number of bins in the region",
          "Bayes Factor",
          "Evidence level (Lee and Wagenmakers, 2013)")
      )
    )

    ## Build Variant data
    vcf<-data.frame()
    if (class(x) == "ExomeDepth") {
      cnv<-x@CNV.calls
      #FILTER
      filter<-as.vector(by(cnv,seq_len(nrow(cnv)), function(c) {
        f<-vector()
        if (abs(c$reads.ratio-1)<=0.25) f<-append(f,"cnvCopyRatio")
        if (c$BF<3) f<-append(f,"cnvLowBF")
        if (is.na(c$exons.hg19)) f<-append(f,"outOfScope")
        if (length(f)==0) f<-append(f,"PASS")
        paste0(f,collapse=',')
      }))
      #REF
      ref<-as.character(
        scanFa(
          fasta,
          GRanges(
            seqnames=cnv$chromosome,
            IRanges(start=cnv$start, end=cnv$start)
          )
        )
      )
      #ALT
      alt<-ifelse(cnv$type=="duplication","<DUP>",ifelse(cnv$type=="deletion","<DEL>","<CNV>"))
      #INFO
      info<-data.frame(
        paste0('END=', cnv$end),
        paste0("SVLEN=", cnv$end-cnv$start),
        paste0("CIPOS=",paste(roi[cnv$start.p-1,"end"]-cnv$start,"0",sep=',')),
        paste0("CIEND=",paste("0",roi[cnv$end.p+1,"start"]-cnv$end,sep=',')),
        paste0('GENES=',sapply(cnv$exons.hg19, function(x) paste(unique(sapply(strsplit(strsplit(ifelse(is.na(x),'.',x),',')[[1]],'_'),function(y) y[1])),collapse=','))),
        paste0("SVTYPE=", ifelse(cnv$type=="duplication","DUP",ifelse(cnv$type=="deletion","DEL","CNV")))
      )
      #FORMAT
      format = data.frame(
        GT=sapply(cnv$reads.ratio, function(x) {
          cn<-round(x*2)
          ifelse(cn>2,'./1',ifelse(cn==2,'./.',ifelse(cn==1,'0/1','1/1')))
        }),
        SM=cnv$reads.ratio,
        CN=round(cnv$reads.ratio*2),
        BC=cnv$end.p-cnv$start.p+1,
        BF=cnv$BF,
        EL=sapply(cnv$BF, function(bf) ifelse(bf>100,'EXTREME',
          ifelse(bf>30,"VERY_STRONG",
            ifelse(bf>10,"STRONG",
              ifelse(bf>3,"MODERATE",
                ifelse(bf>1,"ANECDOTAL","NO_EVIDENCE"))))))
      )
      ## create VCF object
      vcf<-data.frame(
        CHROM=cnv$chromosome,
        POS=cnv$start,
        ID = rep(".",nrow(cnv)), 
        REF = ref,
        ALT = alt,
        QUAL = rep(".",nrow(cnv)),
        FILTER = filter,
        INFO=apply(info,1,function(x) paste(x,collapse=';')),
        FORMAT=rep("GT:SM:CN:BC:BF:EL",nrow(cnv)),
        default=apply(format,1,function(x) paste(x,collapse=':')),
        # samplename
        stringsAsFactors = FALSE
      )
      colnames(vcf)[10]<-samplename
    } else {
      # empty variant dataframe
      vcf<-data.frame(
        CHROM=vector(),
        POS=vector(),
        ID=vector(),
        REF=vector(),
        ALT=vector(),
        QUAL=vector(),
        FILTER=vector(),
        INFO=vector(),
        FORMAT=vector(),
        default=vector(),
        # samplename
        stringsAsFactors = FALSE
      )

    }

    # merge header and cnv data
    vcf <- list(header = header, vcf = vcf)

    # create VCF file
    file.create(filename)

    ## write header
    for (i in 1:length(vcf$header)) {
      if (is.null(dim(vcf$header[[i]]))) {
        # simple
        if (length(vcf$header[[i]]) == 1) {
          # single element
          if (grepl(" ", vcf$header[[i]])) {
            ## has space, add quotes
            vcf$header[[i]] <- dQuote(vcf$header[[i]])
          }
          cat(paste0("##", names(vcf$header)[i], "=", vcf$header[[i]], "\n"), file = filename, append = TRUE)
        } else {
          # muliple elements
          cat(paste0("##", names(vcf$header)[i], "=<", paste(paste(names(vcf$header[[i]]), dQuote(vcf$header[[i]]), sep = "="), collapse = ","), 
                  ">\n"), file = filename, append = TRUE)
        }
      } else {
        # complex (dataframe)
        vcf$header[[i]] <- apply(vcf$header[[i]], 1, function(vcf, name) {
          vcf[grepl(" ", vcf)] <- dQuote(vcf[grepl(" ", vcf)])
          paste0(name, "=", vcf)
        }, colnames(vcf$header[[i]]))
        vcf$header[[i]] <- apply(vcf$header[[i]], 2, function(vcf, name) {
          paste0("##", name, "=<", paste(vcf, collapse = ","), ">")
        }, names(vcf$header)[i])
        # write
        write.table(vcf$header[[i]], filename, quote = FALSE, 
          sep = "", row.names = FALSE, col.names = FALSE, append = TRUE)
        }
    }
    ## write variants with column names
    cat(paste0("#", paste(names(vcf$vcf), collapse = "\t"), "\n"), 
        file = filename, append = TRUE)
    if (!is.null(vcf$vcf)) {
      write.table(vcf$vcf, filename, row.names = FALSE, col.names = FALSE, 
          sep = "\t", quote = FALSE, append = TRUE)
    }
    ## display file for debugging
    # writeLines(readLines(filename))
}
