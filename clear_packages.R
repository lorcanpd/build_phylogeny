#!/usr/bin/env Rscript


installed_packages_info <- installed.packages()[, c("Package", "Version")]
write.csv(installed_packages_info, file = "installed_packages_backup.csv", row.names = FALSE)


# Get a list of all installed packages
installed_packages <- installed.packages()

# Base R packages (as of R 4.1.3, adjust based on your version)
base_packages <- c(
    "base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods",
    "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils"
)

# Filter out base packages from the list of installed packages
additional_packages <- setdiff(installed_packages[, "Package"], base_packages)

# Remove additional packages
remove.packages(additional_packages)
