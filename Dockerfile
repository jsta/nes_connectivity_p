FROM rocker/geospatial:latest
ARG DEBIAN_FRONTEND=noninteractive
ARG LANG=C

MAINTAINER "Joseph Stachelek" stachel2@msu.edu

# RUN dpkg --purge --force-depends tex-common texinfo texlive-local \ 
#   && apt --fix-broken install

# RUN apt-get update \
#  && apt-get install -y --no-install-recommends \
#    apt-utils

# RUN apt-get update \
#  && apt-get install -y --no-install-recommends \
#    texlive-base

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    grass-dev \
    p7zip-full \
    curl \
    libmagick++-dev

RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

RUN tlmgr install beamer translator beamerposter type1cm fp pgfplots subfiles pbox multirow colortbl pdfcrop

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('LAGOSNE','brms','dplyr','tidyr','modelr','data.tree','ggplot2','tidybayes','cowplot','nesRdata','pinp'))"

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_remotes.R');install_remotes(c('jsta/nhdR','jsta/spnetwork','jsta/streamnet','jsta/rgrass7sf'))"

RUN mkdir /liftrroot/
WORKDIR /liftrroot/
