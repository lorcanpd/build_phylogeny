# r_phylogeny
An updated version of the Behjati lab's R pipeline for creating maximum parsimony phylogenies.


## Dependencies
The following R packages need to be installed using bioconductor after 
installing this package:

```{R}
BiocManager::install("GenomicRanges")
BiocManager::install("ggtree")
BiocManager::install("Rsamtools")
```



The following exact dependencies were used to develop the pipeline with `R` version 4.1.3:
- `ape` 5.6-2
- `optparse` 1.7.3
- `seqinr` 4.2-16
- `VGAM` 1.1-7
- `BiocManager` 1.30.18
- `data.table` 1.14.2
- `dplyr` 1.0.9
- `extraDistr` 1.10.0
- `ggplot2` 3.3.6
- `gridExtra` 2.3
- `MASS` 7.3-55
- `doParallel` 1.0.16
- `readr` 2.1.4
- `stringr` 1.4.0
- `tidyr` 1.2.0

To download and install with the exact versions use the following commands:
```
git clone https://github.com/lorcanpd/r_phylogeny.git
Rscript r_phylogeny/install_dependencies.R
```


