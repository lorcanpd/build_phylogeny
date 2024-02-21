# Check that R version 4.1.3
if (as.character(getRversion()) != "4.1.3") {
    stop("R version 4.1.3 is required")
}

# set CRAN mirror to London
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Install specific versions of CRAN packages
install_version <- function(package, version) {
    # Check if devtools is installed, if not install it
    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools", dependencies = TRUE, verbose = TRUE)
    }

    # Use devtools to install specific version
    devtools::install_version(package, version)
}

# Function to compare installed version with desired version
check_version <- function(package, version) {
    installed_version <- as.character(packageVersion(package))
    message(
        sprintf(
            "Package: %s, Installed Version: %s, Expected Version: %s, Match: %s",
            package, installed_version, version, installed_version == version)
    )
}

package_versions <- c(
    ape = "5.6-2",
    optparse = "1.7.3",
    seqinr = "4.2-16",
    VGAM = "1.1-7",
    BiocManager = "1.30.18",
    data.table = "1.14.2",
    dplyr = "1.0.9",
    extraDistr = "1.10.0",
    ggplot2 = "3.3.6",
    gridExtra = "2.3",
    MASS = "7.3-55",
    doParallel = "1.0.16",
    readr = "2.1.4",
    stringr = "1.4.0",
    tidyr = "1.2.0"
)

# Install CRAN packages
lapply(names(package_versions), function(pkg) install_version(pkg, package_versions[pkg]))

# Apply the function to each package
lapply(names(package_versions), function(pkg) check_version(pkg, package_versions[pkg]))


bioconductor_packages <- c(
    GenomicRanges = "1.46.1",
    ggtree = "3.2.1",
    Rsamtools = "2.10.0"
)


lapply(
    names(bioconductor_packages),
    function(pkg) BiocManager::install(
        pkg
    )
)

lapply(
    names(bioconductor_packages),
    function(pkg) check_version(pkg, bioconductor_packages[pkg])
)


