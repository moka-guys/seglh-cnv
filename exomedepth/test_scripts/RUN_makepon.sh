docker run -it -v /home/dbrawand/code/seglh-cnv/data:/data -v /home/dbrawand/code/seglh-cnv/exomedepth:/root -v /srv/work/genome:/genome seglh/exomedepth:latest /root/readCount.R /data/pon/pon.RData /genome/human_g1k_v37_decoy.fasta /data/targets.bed /data/pon/NORMAL_01_XX_M_panel_Pan000_S01_R1_001.bam /data/pon/NORMAL_01_XX_F_panel_Pan000_S01_R1_001.bam /data/pon/NORMAL_02_XX_M_panel_Pan000_S02_R1_001.bam /data/pon/NORMAL_02_XX_F_panel_Pan000_S02_R1_001.bam /data/pon/NORMAL_03_XX_M_panel_Pan000_S03_R1_001.bam /data/pon/NORMAL_03_XX_F_panel_Pan000_S03_R1_001.bam /data/pon/NORMAL_04_XX_M_panel_Pan000_S04_R1_001.bam /data/pon/NORMAL_04_XX_F_panel_Pan000_S04_R1_001.bam /data/pon/NORMAL_05_XX_M_panel_Pan000_S05_R1_001.bam /data/pon/NORMAL_05_XX_F_panel_Pan000_S05_R1_001.bam
