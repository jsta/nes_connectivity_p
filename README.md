[![Paper DOI](https://img.shields.io/badge/Paper-10.1002/lno.11137-blue.svg)](https://doi.org/10.1002/lno.11137) [![Docker Build](https://img.shields.io/docker/build/jsta/stachelek_soranno_2019.svg)](https://cloud.docker.com/repository/docker/jsta/stachelek_soranno_2019)

Code and data for:

**Joseph Stachelek and Patricia A. Soranno (In Press)**. Does Freshwater Connectivity Influence Phosphorus Retention in Lakes? *Limnology and Oceanography* [10.1002/lno.11137](https://doi.org/10.1002/lno.11137)

### Products

Pre-print Manuscript: [manuscript/manuscript.pdf](manuscript/manuscript.pdf)

### Reproducibility

The full system and R environment required to reproduce the paper and analyses is specified in this [Dockerfile](Dockerfile) ([jsta/stachelek_soranno_2019](https://cloud.docker.com/repository/docker/jsta/stachelek_soranno_2019)).

Run the following commands with Docker installed:

```
docker pull jsta/stachelek_soranno_2019
docker run --rm -e PASSWORD=<PASSWORD> jsta/stachelek_soranno_2019
docker ps # note container "code name"
docker exec -ti <NAME> /bin/bash
cd nes_connectivity_p

# build analysis from precomputed connectivity metrics
make -B all

# test connectivity metric calculations
make test_calc_metrics
```

Run the following to recalculate connectivity metrics (requires linking LAGOSNE-GIS as a volume, see https://):

```
docker run -e PASSWORD=<PASSWORD> --rm -v ~/.local/share/LAGOS-GIS:/root/.local/share/LAGOS-GIS jsta/stachelek_soranno_2019
docker ps # note container "code name"
docker exec -ti <NAME> /bin/bash
cd nes_connectivity_p

# test a single lake
Rscript data/connectivity_metrics.R 7003

# recalculate all lakes
make -B data
```
