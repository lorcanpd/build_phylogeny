
# Check that R version 4.1.2
if (R.version$major != 4 | R.version$minor != 1 | R.version$patch != 2) {
    stop("R version 4.1.2 is required")
}

# Install specific versions of CRAN packages
install_version <- function(package, version) {
    # Check if devtools is installed, if not install it
    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
    }

    # Use devtools to install specific version
    devtools::install_version(package, version)
}

# Function to compare installed version with desired version
check_version <- function(package, version) {
    installed_version <- as.character(packageVersion(package))
    message(sprintf("Package: %s, Installed Version: %s, Expected Version: %s, Match: %s",
                    package, installed_version, version, installed_version == version))
}

package_versions <- c(
    tidyr = "1.2.0",
    stringr = "1.4.0",
    data.table = "1.14.2",
    ggplot2 = "3.3.5",
    dplyr = "1.0.7",
    MASS = "7.3-55",
    jsonlite = "1.7.3"
)

# Apply the function to each package
lapply(names(package_versions), function(pkg) check_version(pkg, package_versions[pkg]))

# Install Bioconductor package
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("BiocManager", version = "1.30")

