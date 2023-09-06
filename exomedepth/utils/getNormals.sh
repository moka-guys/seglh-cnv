# Call CNVs from random reference sets (bootstrap) to identify samples that should be excluded from the analysis.
# This example will only assess samples 1 through 10. 
docker run --rm -v `pwd`/data/230824:/home/dnanexus seglh/exomedepth:latest \
	bootstrapCnv.R \
	/home/dnanexus/Pan4398exomedepth_readCount.RData \
	100 30 1 10

# Concatenate score results into a single file.
cat /home/dnanexus/Pan4398exomedepth_readCount.????????????.scores > /home/dnanexus/Pan4398exomedepth_readCount.concatenated.scores

# Display the samples that should be excluded as they show reproducible CNVs with random reference sets.
docker run --rm -v `pwd`/data/230824:/home/dnanexus seglh/exomedepth:latest \
	bootstrapCnv.R \
	/home/dnanexus/Pan4398exomedepth_readCount.RData \
	100 30 \
	/home/dnanexus/Pan4398exomedepth_readCount.concatenated.scores

# rerun the first step with suspected non-normals removed (validation).
docker run --rm -v `pwd`/data/230824:/home/dnanexus seglh/exomedepth:latest \
	bootstrapCnv.R \
	/home/dnanexus/Pan4398exomedepth_readCount.RData \
	100 30 \
	/home/dnanexus/Pan4398exomedepth_readCount.concatenated.scores \
	RUN
