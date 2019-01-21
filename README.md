Code and data for:

**Joseph Stachelek and Patricia A. Soranno (In Press)**. Does Freshwater Connectivity Influence Phosphorus Retention in Lakes? *Limnology and Oceanography*
[preprint](manuscript/manuscript.pdf)

![Docker Build](https://img.shields.io/docker/build/jsta/stachelek_soranno_2019.svg)

### Products

Figures: [manuscript/figures.pdf](manuscript/figures.pdf)

Appendix: [manuscript/appendix.pdf](manuscript/appendix.pdf)

Manuscript: [manuscript/manuscript.pdf](manuscript/manuscript.pdf)

### Reproducibility

The full system and R environment required to reproduce the paper and analyses is specified in this [Dockerfile](Dockerfile) which is also on Dockerhub at [jsta/stachelek_soranno_2019](https://cloud.docker.com/repository/docker/jsta/stachelek_soranno_2019).

Run the following commands with Git and Docker installed:

```
git pull https://github.com/jsta/nes_connectivity_p.git
cd nes_connectivity_p
docker run -e PASSWORD=<PASSWORD> --rm --user root -p 8787:8787 -v /home/jose/test/nes_connectivity_p/:/home/rstudio/stachelek_soranno_2019 jsta/stachelek_soranno_2019
docker ps
docker exec -ti -u root <NAME> /bin/bash
setwd("nes_connectivity_p")
make -B all
```

#### Data requirements

  * LAGOS-NE-GIS
  * `nesRdata`

#### System requirements

* GRASS GIS
  * v.stream.order
  * v.net

* Stan

* R
  * brms
  * LAGOSNE
  * dplyr
  * brms
  * tidybayes
  * ggplot2
  * cowplot
  * data.tree
  * partykit
  * modelr
  * tidyr
  * ggridges
  * nhdR
  * streamnet
  * sf