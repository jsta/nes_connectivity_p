FROM rocker/geospatial:latest
ARG DEBIAN_FRONTEND=noninteractive
ARG LANG=C

MAINTAINER "Jemma Stachelek" stachel2@msu.edu

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    grass-dev \
    p7zip-full \
    curl \
    libmagick++-dev \
    imagemagick \
    pdftk \
    gnupg2

RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

RUN tlmgr install beamer translator beamerposter type1cm fp pgfplots subfiles pbox multirow colortbl pdfcrop

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('LAGOSNE/1.2.0','brms','dplyr','tidyr','modelr','data.tree','ggplot2','cowplot','nesRdata','pinp','pheatmap','USAboundaries','partykit','ggsn','vapour','kableExtra'))"

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_remotes.R');install_remotes(c('jsta/nhdR','jsta/spnetwork','jsta/streamnet','jsta/rgrass7sf','jsta/tidybayes','drsimonj/corrr@v0.2.1','jsta/LAGOSNEgis','jsta/nlaR'))"

RUN Rscript -e "dir.create(rappdirs::user_data_dir(), recursive = TRUE)"

RUN Rscript -e "streamnet:::install_grass_extensions()"

RUN Rscript -e "LAGOSNE::lagosne_get('1.087.1', dest_folder = LAGOSNE:::lagos_path())"

RUN Rscript -e "nlaR::nla_get(2012, use_rappdirs = TRUE)"

RUN git clone https://github.com/jsta/nes_connectivity_p.git
