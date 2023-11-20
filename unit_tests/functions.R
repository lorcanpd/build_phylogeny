require(testthat)

# Tests for plot_spectrum function
source("../functions.R")

# TODO replace placeholder test with real tests.
# replace the expected values with the actual values expected based on the mock
# data.

# Mock data for testing
mock_data <- data.frame(
    chr = c("1", "2"),
    pos = c(100, 200),
    ref = c("A", "C"),
    mut = c("G", "T")
)

# Test data processing
test_that("Data is processed correctly", {
    processed_data <- plot_spectrum(
        bed = mock_data, save = NULL, genomeFile = "path/to/genome.fa"
    )
    # Replace with actual expected values
    expect_equal(nrow(processed_data$mutations), 2)
    expect_true(all(processed_data$mutations$ref %in% c("A", "C", "G", "T")))
    # Add more expectations as needed
})

# Test frequency calculation
test_that("Frequency calculation is correct", {
    freq_data <- plot_spectrum(
        bed = mock_data, save = NULL, genomeFile = "path/to/genome.fa"
    )
    # Replace with actual expected values
    expect_equal(length(freq_data$freqs_full), expected_length)
    # Add more expectations as needed
})

