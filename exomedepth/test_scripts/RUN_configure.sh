docker run -it \
	-v /home/dbrawand/code/seglh-cnv/data:/data \
	-v /home/dbrawand/code/seglh-cnv/exomedepth:/test \
	-v /srv/work/genome:/genome \
	seglh/exomedepth:latest \
	/test/configure.R \
	/data/pon/qc.RData \
	/data/pon/annot.bed
