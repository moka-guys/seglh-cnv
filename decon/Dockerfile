# Base Image
FROM r-base:3.1.2

# Metadata
LABEL base.image="r-base:r3.1.2"
LABEL software="DeCON"
LABEL description="CNV calling for DNA capture assays"

# Maintainer
MAINTAINER David Brawand <dbrawand@nhs.net>

# set working directory
WORKDIR /root
ENV HOME /root

RUN apt-get -y update
RUN apt-get --yes install libcurl4-openssl-dev
RUN apt-get --yes install time
RUN apt-get --yes install git

# install DECoN
RUN git clone --recursive https://github.com/RahmanTeam/DECoN.git
WORKDIR /root/DECoN/Linux
RUN sed -e '/options/ s/^#*/#/' -i packrat/packrat_source/.Rprofile
RUN ./setup.sh
RUN Rscript sessionInfo.R --bootstrap-packrat

# LateX packages
#RUN tlmgr init-usertree && \
#	tlmgr update --self && \
#	tlmgr install framed
#
## R dependencies
#RUN Rscript -e "install.packages('BiocManager');BiocManager::install(version = '3.12')" && \
#	Rscript -e "BiocManager::install(c('Biostrings','IRanges','Rsamtools','GenomicRanges','GenomicAlignments'))" && \
#	Rscript -e "install.packages('optparse')" && \
#	Rscript -e "install.packages('ExomeDepth')" && \
#	Rscript -e "install.packages('plyr')" && \
#	Rscript -e "install.packages('knitr')" && \
#	Rscript -e "install.packages('tinytex')"
#
## add exomedepth scripts and report template
#ADD readCount.R /root
#ADD exomeDepth.R /root
#ADD exomeDepth.Rnw /root

ENTRYPOINT ["/usr/bin/Rscript"]
