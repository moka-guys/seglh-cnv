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

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    time \
    libpoppler-cpp-dev \
    libssl-dev \
    libfontconfig1-dev \
    libxml2-dev \
    wget \
    perl \
    xz-utils \
    && rm -rf /var/lib/apt/lists/*

# Copy TeX Live profile
COPY texlive.profile /tmp/texlive.profile

# Install TeX Live 2020
RUN wget https://pi.kwarc.info/historic/systems/texlive/2020/install-tl-unx.tar.gz && \
    tar -xzf install-tl-unx.tar.gz && \
    cd install-tl-2020* && \
    ./install-tl -profile /tmp/texlive.profile -repository ftp://tug.org/historic/systems/texlive/2020/tlnet-final && \
    cd .. && rm -rf install-tl-unx.tar.gz install-tl-2020*

# Add TeX Live to PATH
ENV PATH=/usr/local/texlive/2020/bin/x86_64-linux:$PATH

# Initialize tlmgr and set repository to 2020
RUN tlmgr init-usertree && \
    tlmgr option repository ftp://tug.org/historic/systems/texlive/2020/tlnet-final && \
    tlmgr update --self
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
	tlmgr install lastpage && \
	tlmgr install pgf && \
	tlmgr install xcolor

# set CRAN mirror
RUN echo 'options(repos=structure(c(CRAN="http://cran.ma.imperial.ac.uk/")))' > /root/.Rprofile && \
    echo 'options(BioC_mirror="http://bioconductor.statistik.tu-dortmund.de/")' >> /root/.Rprofile

# R packages (BioC)
RUN Rscript -e "install.packages('BiocManager');BiocManager::install(version = '3.15', ask = FALSE)"
RUN Rscript -e "BiocManager::install(c('Biostrings','IRanges','Rsamtools','GenomicRanges','GenomicAlignments'))"

# R packages (CRAN)
RUN Rscript -e "install.packages('optparse')" && \
	Rscript -e "install.packages('randomForest')" && \
	Rscript -e "install.packages('aod')" && \
	Rscript -e "install.packages('ExomeDepth')" && \
	Rscript -e "install.packages('plyr')" && \
	Rscript -e "install.packages('knitr')" && \
	Rscript -e "install.packages('kableExtra')" && \
	Rscript -e "install.packages('tinytex')" && \
	Rscript -e "install.packages('xtable')" && \
	Rscript -e "install.packages('moments')"

# install ED from source as CRAN rejected ED (10/2022) v1.1.16
# CRAN reinstated ED (11/2020) v1.1.16
#RUN Rscript -e "install.packages('devtools')" && \
#    Rscript -e "require(devtools);install_github('vplagnol/ExomeDepth', ref='7f6f8621801d1eb46880f88c4f71bb7c31680c73')"
RUN apt-get update && apt-get install -y \
    libcairo2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    pandoc \
    git \
    make \
    g++ \
    && rm -rf /var/lib/apt/lists/*
# Install kableExtra from github as CRAN package not working
RUN R -e "install.packages(c('remotes','svglite','textshaping'), repos='http://cran.ma.imperial.ac.uk/')"
RUN R -e "remotes::install_github('haozhu233/kableExtra')"
# Install Exomdepth from github as CRAN package not working
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org'); \
                remotes::install_github('vplagnol/ExomeDepth', ref='7f6f8621801d1eb46880f88c4f71bb7c31680c73')"
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
