# Base Image
FROM rstudio/r-base:4.2.1-focal

# Metadata
LABEL base.image="r-base:4.2.1"
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
RUN apt-get --yes install \
	libcurl4-openssl-dev \
	time \
	libpoppler-cpp-dev \
	texlive \
	libssl-dev \
	libfontconfig1-dev \
	libxml2-dev

# LateX packages (Locked on 2020 release)
RUN tlmgr init-usertree
# RUN tlmgr repository add ftp://tug.org/historic/systems/texlive/2020/tlnet-final
# RUN tlmgr repository remove https://mirror.ctan.org/systems/texlive/tlnet
# RUN tlmgr option repository ftp://tug.org/historic/systems/texlive/2020/tlnet-final
# RUN tlmgr update --self
RUN tlmgr install framed && \
	tlmgr install mdframed && \
	tlmgr install zref && \
	tlmgr install needspace && \
	tlmgr install dingbat && \
	tlmgr install booktabs && \
	tlmgr install makecell && \
	tlmgr install colortbl && \
	tlmgr install parskip && \
	tlmgr install fancyhdr && \
	tlmgr install xcolor

# set CRAN mirror
RUN echo 'options(repos=structure(c(CRAN="http://cran.ma.imperial.ac.uk/")))' > /root/.Rprofile

# R packages (BioC)
RUN Rscript -e "install.packages('BiocManager');BiocManager::install(version = '3.15')" && \
	Rscript -e "BiocManager::install(c('Biostrings','IRanges','Rsamtools','GenomicRanges','GenomicAlignments'))"

# R packages (CRAN)
RUN Rscript -e "install.packages('optparse')" && \
	Rscript -e "install.packages('randomForest')" && \
	Rscript -e "install.packages('aod')" && \
	Rscript -e "install.packages('ExomeDepth')" && \
	Rscript -e "install.packages('plyr')" && \
	Rscript -e "install.packages('knitr')" && \
	Rscript -e "install.packages('kableExtra')" && \
	Rscript -e "install.packages('tinytex')" && \
	Rscript -e "install.packages('xtable')"
#RUN Rscript -e "install.packages('Hmisc')"

# add exomedepth scripts and report template
ADD readCount.R /root
ADD ed2vcf.R /root
ADD exomeDepth.R /root
ADD exomeDepth.Rnw /root
ADD VERSION /root

# add QC configs
COPY qc_config/* /root/

# Rscript as entrypoint
ENTRYPOINT ["/usr/bin/Rscript"]
