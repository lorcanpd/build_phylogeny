# build_phylogeny
This is a rewrite of https://github.com/TimCoorens/Sequoia. It has been modularised and rewritten for readability, re-usability, extensibility, and efficiency.


## Use

---
Please see `example.R` for an end-to-end example of how to use the package.
`example_params.json` contains the parameters used by the example.



## Installation

---

### Install using `devtools`

```{R}
devtools::install_github("lorcanpd/r_phylogeny")
devtools::install_git('https://github.com/NickWilliamsSanger/treemut')

# The following R packages need to be installed using bioconductor after 
installing this package:

BiocManager::install("GenomicRanges")
BiocManager::install("ggtree")
BiocManager::install("Rsamtools")
```
Followed by:
```{bash}
# Install mpboot
git clone https://github.com/diepthihoang/mpboot.git && \
    mkdir /mpboot-build && \
    cd /mpboot-build && \
    cmake ../mpboot -DIQTREE_FLAGS=avx -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ && \
    make -j4

# Ensure mpboot executable is in the PATH
ENV PATH="/mpboot-build:$PATH"
```
---
### Cloning the repository

To download and install with the exact dependencies used indevelopment use the following commands:
```{bash}
git clone https://github.com/lorcanpd/r_phylogeny.git
Rscript r_phylogeny/install_dependencies.R
R CMD INSTALL r_phylogeny/.
```



