[![Paper DOI](https://img.shields.io/badge/Paper-10.1002/lno.11137-blue.svg)](https://doi.org/10.1002/lno.11137) [![Code DOI](https://zenodo.org/badge/123951266.svg)](https://zenodo.org/badge/latestdoi/123951266) [![Docker Build](https://img.shields.io/badge/Docker%20Image-jsta/stachelek--soranno--2019-green.svg)](https://cloud.docker.com/repository/docker/jsta/stachelek_soranno_2019)

Code and data for:

**Joseph Stachelek and Patricia A. Soranno. 2019**. Does Freshwater Connectivity Influence Phosphorus Retention in Lakes? *Limnology and Oceanography* [10.1002/lno.11137](https://doi.org/10.1002/lno.11137)

### Reproducibility

The full system and R environment required to reproduce the paper and analyses is specified in this [Dockerfile](Dockerfile) ([jsta/stachelek_soranno_2019](https://cloud.docker.com/repository/docker/jsta/stachelek_soranno_2019)).

Run the following commands with Docker installed to rebuild the pre-print and test connectivity metric calculations:

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

Run the following commands to fully recalculate connectivity metrics:

  * requires linking to LAGOSNE-GIS as a Docker volume, get data with the `LAGOSNEgis` package [![DOI](https://zenodo.org/badge/106293356.svg)](https://zenodo.org/badge/latestdoi/106293356)

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

#### Bespoke (non-CRAN) dependencies

|||
|---|---|
|`LAGOSNEgis`| [![DOI](https://zenodo.org/badge/106293356.svg)](https://zenodo.org/badge/latestdoi/106293356)|
|`nhdR`| [![DOI](https://zenodo.org/badge/75339263.svg)](https://zenodo.org/badge/latestdoi/75339263)|
|`spnetwork`| [![DOI](https://zenodo.org/badge/96037556.svg)](https://zenodo.org/badge/latestdoi/96037556)|
|`streamnet`| [![DOI](https://zenodo.org/badge/104792308.svg)](https://zenodo.org/badge/latestdoi/104792308)|
|`rgrass7sf`| [![DOI](https://zenodo.org/badge/115946587.svg)](https://zenodo.org/badge/latestdoi/115946587)|
|`tidybayes`| [![DOI](https://zenodo.org/badge/116701609.svg)](https://zenodo.org/badge/latestdoi/116701609)|
|`nlaR`| [![DOI](https://zenodo.org/badge/75324775.svg)](https://zenodo.org/badge/latestdoi/75324775)|
