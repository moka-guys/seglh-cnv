docker run -it \
	-v /srv/work/genome:/genome \
	-v /home/dbrawand/code/seglh-cnv/data:/data \
	seglh/decon \
	ReadInBams.R \
	--bams /data \
	--bed /data/targets.decon.bed \
	--fasta /genome/human_g1k_v37_decoy.fasta \
	--out /data/decon_counts

docker run -it \
	-v /srv/work/genome:/genome \
	-v /home/dbrawand/code/seglh-cnv/data:/data \
	seglh/decon \
	IdentifyFailures.R \ 
	--Rdata /data/decon_counts.RData \
	--mincorr 0.98 \
	--mincov 100 \
	--output /data/decon_qc

