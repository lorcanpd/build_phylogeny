# build_phylogeny
This is a rewrite of https://github.com/TimCoorens/Sequoia. It has been modularised and rewritten for readability, re-usability, and extensibility. Functions have been rewritten to take advantage of parallel processing for greater efficiency.

---

## Use


Please see `example.R` for an end-to-end example of the package being used to construct a phylogeny. `example_params.json` contains the parameters used by the example. The default  parameters (`default_params`) are described inside the `read_params.R` module.

---

## Installation



### Create `Docker` image and convert to `Singularity` image


Copy the `Dockerfile.dev` from this repository to a new directory and build the image using the following command:

```{bash}
docker build --tag build_phylogeny --file Dockerfile.dev .
```
Then convert the `Docker` image to a `Singularity` image using the following command:

```{bash}
singularity build build_phylogeny.sif docker-daemon://build_phylogeny:latest
```
Once the image is created it can be used to run the `build_phylogeny` package.

---

### Install using `devtools`

If you'd rather not make a container... 

```{R}
devtools::install_git("https://github.com/NickWilliamsSanger/treemut")

# The following R packages need to be installed using bioconductor after 
installing this package:

BiocManager::install("GenomicRanges")
BiocManager::install("ggtree")
BiocManager::install("Rsamtools")

devtools::install_github("https://github.com/lorcanpd/build_phylogeny")
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




