#!/usr/bin/sh

# command lines to run exomedepth analysis
#
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


# X68513f74:NG11-131-5:0.0001
# +-key      |         +-transition probability threshold (optional)
#            +-human readable name
#
docker run -it \
	-v /path_to_data:/data \
	seglh/exomedepth:latest \
	exomeDepth.R \
	VERSIONSTRING \
	/data/output.pdf \
	/data/targets.bed:PANELNAME \
	/data/readCount.RData \
	/data/sample1.bam:Sample1

