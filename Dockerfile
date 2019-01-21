FROM rocker/geospatial:latest

MAINTAINER Joseph Stachelek <stachel2@msu.edu>

# System dependencies for required R packages
RUN  rm -rf  /var/cache/apt/* \
  && apt-get update -qq \
  && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    ca-certificates \
    apt-utils \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    git

# RUN dpkg --purge --force-depends tex-common && apt -y --fix-broken install

RUN apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends grass-dev p7zip-full curl texlive-base texlive-extra-utils 

RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

RUN tlmgr install beamer translator beamerposter

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('LAGOSNE','brms','dplyr','tidyr','modelr','data.tree','ggplot2','tidybayes','cowplot','nesRdata', 'pinp'))"

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_remotes.R');install_remotes(c('jsta/nhdR','jsta/spnetwork','jsta/streamnet','jsta/rgrass7sf'))"

RUN mkdir /liftrroot/
WORKDIR /liftrroot/
