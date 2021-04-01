---
title: "ExomeDepth CNV pipeline"
author: "David Brawand <dbrawand@nhs.net>"
date: "2021/02/15"
output:
  pdf_document:
    toc: true
    highlight: zenburn
---

# SEGLH ExomeDepth pipeline

## input files

Exomdepth runs on a collection of intervals that must be provided in BED format.
The target intervals must be annotated with gene symbol and exon number as follows:

```
3	128199861	128200161	GATA2_5
3	128200661	128200787	GATA2_4
3	128202702	128202848	GATA2_3
3	128204569	128205211	GATA2_2
3	128205645	128205874	GATA2_1
```
If intronic or intergenic regions are included they can be named accordingly using the same naming principle (eg. GATA2-INTRON_1). Exomedepth uses the exon/intron numbers only to label the CNVs in the report output and they are not part of the applied mathematical model.

## Usage

Build the docker container which contains all scripts with `make`.
The built docker image is tagged as `seglh/exomedepth:latest`.

The docker image contains all required packages and scripts but no reference genome.

### Step 1 - Create read count data

```
docker run -it \
	-v /path_to_data:/data \
	-v /path_to_genome:/resources \
	seglh/exomedepth:latest \
	readCount.R \
	/data/readCount.RData \
	/genome/human_g1k_v37_decoy.fasta \
	/data/targets.bed \
	/data/sample1.bam \
	/data/sample2.bam \
	/data/sample3.bam \
	/data/sample4.bam
```

### Step 2 - Run Exomedepth

This steps will pick the most appropriate reference samples and run exomedepth for a given sample.
Reference samples are picked from within the same batch. If any BAM file is prefixed with _NNN_ it will be considered part of a panel of normals which will be used as reference instead of other batched samples (requires at least 2 normal samples).

```
docker run -it \
	-v /path_to_data:/data \
	seglh/exomedepth:latest \
	exomeDepth.R \
	VERSIONSTRING \
	/data/output.pdf \
	/data/targets.bed:PANELNAME \
	/data/readCount.RData \
	/data/sample1.bam:SampleLabel \
	/data/qcmodel.RData
```
NB: The last argument is optional (see below)

#### Output files

1. `output.pdf` contains the CNV report
2. `output.bed` contains the deletion/duplication intervals and the log-likelihood ratio as score.

### Integrated Quality Control

The report contains pertinent information for Quality control. Two plots visualise the RPKM correlation with other samples in the same batch, and the coefficient of variation.

A table summarises if QC thresholds have been met (PASS), failed a soft threshold (CAUTION) or a hard threshold (FAIL). If a threshold is not met, an interpretation help is provided in the report.

Optionally a random forest classifier can be provided which will indicate the PASS/FAIL classification in the report. This classification is only a recommendation and the scientist has ultimate authority to classify the sample as failed.
#### QC Thresholds

Exomdedepth implements default QC thresholds as follows:

| Metric                              | Soft Threshold | Hard Threshold |
| :---------------------------------- | -------------: | -------------: |
| Median Correlation in Batch         |             NA |           0.90 |
| Maximum Correlation in Batch        |           0.95 |           0.90 |
| Coefficent of Variation             |             30 |             35 |
| Correlation with selected reference |           0.95 |           0.90 |
| Selected reference samples          |              5 |              3 |

These can be overridden by providing an RData file containing thresholds as follows:

```
limits<-list(
  medcor=c(NA, 0.90),    # median correlation within batch
  maxcor=c(0.95, 0.90),  # max correlation within batch
  refcor=c(0.95, 0.90),  # reference set correlation
  refcount=c(5, 3),      # refernce set size (selected reference samples)
  coeffvar=c(30, 35)     # coefficient of variation
)
```

The `buildQcRfc.R` script contains an example how these can be specified.

#### Random Forest QC classifier

SEGLH-Exomedepth can optionally provide a QC classification based on previous labeled data. The classfier is automatically optimised for the number of random variables considered and runs 200 trees by default. All labeled samples are part of the training set.


The model can be built as follows:

1. Amend a column with header `status` to the CSV metrics output from the `readCount.R` step.
2. Label the samples with a known classification (eg. FAIL or PASS) in this column. Providing a label is optional. Missing values will not be included in the model generation.
3. Build the model with `Rscript buildQcRfc.R qcmodel.RData file1.csv [file2.csv ...]`
4. Supply the `qcmodel.RData` file to the exomedepth command line as last parameter. The report will automaticcly be amended with the RF QC classification and an explanation.

It is recommended to inspect the generated Random Forest classifier in terms of error rate convergence and out-of-bag error. Also consider adjusting of the number of trees!
