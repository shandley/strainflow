pkgs <- c("testthat", "yaml", "readr", "dplyr", "tibble", "glue", "fs",
           "rlang", "stringr", "processx", "jsonlite", "tidyr", "purrr",
           "ggplot2", "scales", "viridis", "ggrepel", "survival", "lme4",
           "broom", "patchwork", "cli")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cran.r-project.org")
  }
}
cat("All required packages are available.\n")
