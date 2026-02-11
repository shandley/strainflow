# Source all R files and run tests
library(testthat)

# Source all package R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in r_files) {
  tryCatch(
    source(f, local = TRUE),
    error = function(e) message("Warning sourcing ", f, ": ", e$message)
  )
}

# Run tests
test_results <- test_dir("tests/testthat", reporter = "summary")
