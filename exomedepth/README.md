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
``````

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
	/data/sample1.bam:Sample1
``````

#### Output files

1. `output.pdf` contains the CNV report
2. `output.bed` contains the deletion/duplication intervals and the log-likelihood ratio as score.


