docker run -it \
	-v /home/dbrawand/code/seglh-cnv/data:/data \
	-v /home/dbrawand/code/seglh-cnv/exomedepth:/test \
	-v /srv/work/genome:/genome \
	seglh/exomedepth:latest \
	/test/exomeDepth.R \
	DEV \
	/data/pon/output.pdf \
	/data/targets.bed:FULL \
	/data/pon/counts.RData \
	/data/pon/NGS000_00_XX_M_panel_Pan000_S00_R1_001.bam:MALE
