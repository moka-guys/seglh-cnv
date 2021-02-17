# Base Image
FROM r-base:4.0.3

# Metadata
LABEL base.image="r-base:4.0.3"
LABEL software="ExomeDepth"
LABEL description="CNV calling for DNA capture assays"
LABEL website="https://CRAN.R-project.org/package=ExomeDepth"
LABEL documentation="https://CRAN.R-project.org/package=ExomeDepth"
LABEL license="https://CRAN.R-project.org/package=ExomeDepth"

# Maintainer
MAINTAINER David Brawand <dbrawand@nhs.net>

# set working directory
WORKDIR /root
ENV HOME /root

RUN apt-get -y update
RUN apt-get --yes install libcurl4-openssl-dev
RUN apt-get --yes install time
RUN apt-get --yes install libpoppler-cpp-dev
RUN apt-get --yes install texlive

# LateX packages
RUN tlmgr init-usertree && \
	tlmgr update --self && \
	tlmgr install framed

# R dependencies
RUN Rscript -e "install.packages('BiocManager');BiocManager::install(version = '3.12')" && \
	Rscript -e "BiocManager::install(c('Biostrings','IRanges','Rsamtools','GenomicRanges','GenomicAlignments'))" && \
	Rscript -e "install.packages('optparse')" && \
	Rscript -e "install.packages('ExomeDepth')" && \
	Rscript -e "install.packages('plyr')" && \
	Rscript -e "install.packages('knitr')" && \
	Rscript -e "install.packages('tinytex')"

# add exomedepth scripts and report template
ADD readCount.R /root
ADD exomeDepth.R /root
ADD exomeDepth.Rnw /root

ENTRYPOINT ["/usr/bin/Rscript"]
