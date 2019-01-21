FROM rocker/geospatial:latest
ARG DEBIAN_FRONTEND=noninteractive

MAINTAINER "Joseph Stachelek" stachel2@msu.edu

RUN dpkg --purge --force-depends tex-common texinfo \ 
  && apt --fix-broken install

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    apt-utils

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    texlive-base

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    grass-dev \
    p7zip-full \
    curl \
    texlive-latex-extra

RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

RUN tlmgr install beamer translator beamerposter

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('LAGOSNE','brms','dplyr','tidyr','modelr','data.tree','ggplot2','tidybayes','cowplot','nesRdata','pinp'))"

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_remotes.R');install_remotes(c('jsta/nhdR','jsta/spnetwork','jsta/streamnet','jsta/rgrass7sf'))"

RUN mkdir /liftrroot/
WORKDIR /liftrroot/
