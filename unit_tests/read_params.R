library(testthat)
source("../read_params.R")

# Define a test JSON file path (you will need to create this file for testing)
test_json_file <- "test_params.json"

# TODO: fix this up.

# Test 1: Function correctly reads and merges parameters from JSON file
test_that("read and merge parameters from JSON", {
    # Create a test JSON file with some parameters
    test_params <- list(donor_id = "TestPatient", ncores = 4)
    write_json(test_params, test_json_file)

    # Read parameters using your function
    params <- read_json_params(test_json_file, default_parameters)

    # Check if parameters are read and merged correctly
    expect_equal(params$donor_id, "TestPatient")
    expect_equal(params$ncores, 4)
    expect_equal(params$min_cov, default_parameters$min_cov) # Unchanged default parameter
})

# Test 2: Function handles missing JSON file
test_that("handle missing JSON file", {
    # Attempt to read parameters from a non-existent JSON file
    params <- read_json_params("nonexistent.json", default_parameters)

    # Check if default parameters are returned
    expect_equal(params, default_parameters)
})

# Test 3: Function handles extra parameters in JSON file
test_that("handle extra parameters in JSON", {
    # Create a test JSON file with an extra parameter
    test_params <- list(extra_param = "extra")
    write_json(test_params, test_json_file)

    # Read parameters using your function
    params <- read_json_params(test_json_file, default_parameters)

    # Check if extra parameter is ignored
    expect_false("extra_param" %in% names(params))
})

# Clean up: Remove test JSON file
unlink(test_json_file)
