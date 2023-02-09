###########################################################################
## Copyright 2018 Nantomics LLC                                          ##
##                                                                       ##
## Redistribution and use in source and binary forms, with or without    ##
## modification, are permitted for educational, research and non-profit  ##
## purposes, by non-profit institutions only provided that the following ##
## conditions are met:                                                   ##
## 1. Redistributions of source code must retain the above copyright     ##
## notice, this list of conditions and the following disclaimer.         ##
##                                                                       ##
## 2. Redistributions in binary form must reproduce the above copyright  ##
## notice, this list of conditions and the following disclaimer in the   ##
## documentation and/or other materials provided with the distribution.  ##
##                                                                       ##
## 3. Neither the name of the copyright holder nor the names of its      ##
## contributors may be used to endorse or promote products derived from  ##
## this software without specific prior written permission.              ##
##                                                                       ##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS   ##
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ##
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ##
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT  ##
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,            ##
## INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,  ##
## BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS ##
## OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED    ##
## AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT           ##
## LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY ##
## WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE           ##
## POSSIBILITY OF SUCH DAMAGE.                                           ##
###########################################################################

FROM r-base:latest

MAINTAINER Winston Chang "winston@rstudio.com"

# Install dependencies and Download and install shiny server
RUN apt-get update && apt-get install -y -t unstable \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    openjdk-8-jdk \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev && \
    wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb

RUN Rscript -e "install.packages(c('shiny', 'rmarkdown', 'devtools', 'ggplot2', 'dplyr', 'ggvis', 'feather'), repos='https://cran.rstudio.com/')" && \
    Rscript -e "library(devtools); install_github('igraph/rigraph'); install_github('satijalab/seurat')" && \
    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
    rm -rf /var/lib/apt/lists/*

EXPOSE 3838

COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]