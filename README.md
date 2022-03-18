## Research compendium for ‘Carbon and oxygen stable isotopic evidence for diverse sheep and goat husbandry strategies amid a Final Bronze Age farming milieu in the Kyrgyz Tian Shan’

### Compendium DOI:

<http://dx.doi.org/10.17605/OSF.IO/UTNYR>

The files at the URL above will generate the results as found in the publication. The files hosted at <https://github.com/nevrome/isotopes.sheep.goat.tianshan> are the development versions and may have changed since the paper was published.

### Authors of this repository:

- Taylor Hermes [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--8377--468X-green.svg)](http://orcid.org/0000-0002-8377-468X)
- Clemens Schmid [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0003--3448--5715-green.svg)](http://orcid.org/0000-0003-3448-5715)

### Published in:

**(in Print)**: *International Journal of Osteoarchaeology*

### Overview of contents:

This repository contains the directories `code` and `data` to reproduce the data preparation, calculations and figure renderings in this paper. The `plots` and `tables` directories contain readily rendered versions of plots and tables.

`fetch_data.R` downloads, transforms and collates raw input data (`data/raw/...`) for this and from previous publications into a useful and unified table format (`data/input/...`). As we stores these technically intermediate data products, it is generally not necessary to run `fetch_data.R` to reproduce the analysis in this paper.

The actual computation happens in `main_workflow.R`, from which we outsourced various bigger code chunks into separate files. They are sourced automatically by the script when needed. `map.R` is independent of this workflow and only renders the map in Figure 1.

### How to reproduce:

As the data and code in this repository are complete and self-contained, it can be reproduced with only an R environment (tested for R v4.1.0). The necessary package dependencies are documented in the `DESCRIPTION` file and can be installed manually or automatically with 

```r
devtools::install(repos = "https://mran.microsoft.com/snapshot/2022-03-18")
```

### Licenses:

[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/) year: 2022, copyright holder: Taylor Hermes
