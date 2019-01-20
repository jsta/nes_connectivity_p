Code and data for:

**Joseph Stachelek and Patricia A. Soranno (In Press)**. Does Freshwater Connectivity Influence Phosphorus Retention in Lakes? *Limnology and Oceanography*
[preprint](manuscript/manuscript.pdf)

> Limnology and Oceanography

### Products

Figures: [manuscript/figures.pdf](manuscript/figures.pdf)

Appendix: [manuscript/appendix.pdf](manuscript/appendix.pdf)

Manuscript: [manuscript/manuscript.pdf](manuscript/manuscript.pdf)

### Reproducibility

The full system and R environment required to reproduce the paper and analyses is specified in [Dockerfile](Dockerfile) and on Dockerhub [jsta/stachelek_soranno_2019](https://cloud.docker.com/repository/docker/jsta/stachelek_soranno_2019).

Run `docker run -e PASSWORD=<PASSWORD> --rm -p 8787:8787 jsta/stachelek_soranno_2019`
Run `make -B all` to reproduce outputs.

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