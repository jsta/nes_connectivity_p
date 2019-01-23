[![Paper DOI](https://img.shields.io/badge/Paper-10.1002/lno.11137-blue.svg)](https://doi.org/10.1002/lno.11137) [![Docker Build](https://img.shields.io/docker/build/jsta/stachelek_soranno_2019.svg)](https://cloud.docker.com/repository/docker/jsta/stachelek_soranno_2019)

Code and data for:

**Joseph Stachelek and Patricia A. Soranno (In Press)**. Does Freshwater Connectivity Influence Phosphorus Retention in Lakes? *Limnology and Oceanography* [10.1002/lno.11137](https://doi.org/10.1002/lno.11137)

### Products

Figures: [manuscript/figures.pdf](manuscript/figures.pdf)

Appendix: [manuscript/appendix.pdf](manuscript/appendix.pdf)

Pre-print Manuscript: [manuscript/manuscript.pdf](manuscript/manuscript.pdf)

### Reproducibility

The full system and R environment required to reproduce the paper and analyses is specified in this [Dockerfile](Dockerfile) ([jsta/stachelek_soranno_2019](https://cloud.docker.com/repository/docker/jsta/stachelek_soranno_2019)).

Run the following commands with Git and Docker installed:

```
git pull https://github.com/jsta/nes_connectivity_p.git
cd nes_connectivity_p
docker pull jsta/stachelek_soranno_2019
docker run -e PASSWORD=<PASSWORD> --rm -v /home/jose/test/nes_connectivity_p/:/home/jose/nes_connectivity_p -v ~/.local/share/LAGOS-GIS:/home/jose/.local/share/LAGOS-GIS -v ~/Documents/Science/Data/lagos_gis/:/home/jose/Documents/Science/Data/lagos_gis/ jsta/stachelek_soranno_2019
docker ps
docker exec -ti --user root <NAME> /bin/bash
setwd("nes_connectivity_p")

# test metrics calculations
make test

# build all
make -B all
```

