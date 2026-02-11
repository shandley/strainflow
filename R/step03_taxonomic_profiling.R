# ---- step03_taxonomic_profiling.R ----
# Taxonomic profiling with MetaPhlAn4 (Python tool, runs via conda).

#' Run MetaPhlAn4 taxonomic profiling
#'
#' Executes MetaPhlAn4 inside a conda environment to produce a taxonomic
#' profile from non-host metagenomic reads.
#'
#' @param sample_id Character scalar. Sample identifier.
#' @param nonhost_r1 Character scalar. Path to non-host forward reads.
#' @param nonhost_r2 Character scalar. Path to non-host reverse reads.
#' @param output_dir Character scalar. Output directory (created if needed).
#' @param config Named list. Pipeline configuration. Expected elements under
#'   `config$metaphlan`:
#'   \describe{
#'     \item{conda_env}{Character. Name of the conda environment with
#'       MetaPhlAn4 installed (default `"metaphlan4"`).}
#'     \item{threads}{Integer. Number of processing threads (default 8).}
#'     \item{analysis_type}{Character. MetaPhlAn analysis type
#'       (default `"rel_ab_w_read_stats"`).}
#'     \item{database}{Character. Path to MetaPhlAn database (optional).
#'       If provided, `--bowtie2db` is appended to the command.}
#'   }
#'
#' @return A single-row [tibble::tibble] with columns: `sample_id`,
#'   `profile_path`, `bowtie2_out_path`.
#'
#' @export
run_metaphlan4 <- function(sample_id, nonhost_r1, nonhost_r2,
                           output_dir, config) {

  # -- Validate inputs --------------------------------------------------------
  if (!fs::file_exists(nonhost_r1)) {
    rlang::abort(glue::glue("Non-host R1 not found: {nonhost_r1}"))
  }
  if (!fs::file_exists(nonhost_r2)) {
    rlang::abort(glue::glue("Non-host R2 not found: {nonhost_r2}"))
  }
  fs::dir_create(output_dir)

  # -- Resolve config ---------------------------------------------------------
  mp_cfg <- config$metaphlan %||% list()
  conda_env     <- mp_cfg$conda_env     %||% "metaphlan4"
  threads       <- mp_cfg$threads        %||% 8L
  analysis_type <- mp_cfg$analysis_type  %||% "rel_ab_w_read_stats"
  database      <- mp_cfg$database       %||% NULL

  # -- Build output paths -----------------------------------------------------
  profile_path    <- fs::path(output_dir, glue::glue("{sample_id}_metaphlan4.tsv"))
  bowtie2_out     <- fs::path(output_dir, glue::glue("{sample_id}_metaphlan4_bt2.bz2"))

  # -- Build MetaPhlAn command ------------------------------------------------
  input_files <- glue::glue("{nonhost_r1},{nonhost_r2}")

  cmd_parts <- glue::glue(
    "metaphlan {input_files} ",
    "--input_type fastq ",
    "--output_file {profile_path} ",
    "--bowtie2out {bowtie2_out} ",
    "--nproc {threads} ",
    "-t {analysis_type} ",
    "--sample_id {sample_id}"
  )

  # Optionally specify database path

  if (!is.null(database) && nzchar(database)) {
    cmd_parts <- glue::glue("{cmd_parts} --bowtie2db {database}")
  }

  run_conda_cmd(
    conda_env = conda_env,
    cmd       = cmd_parts
  )

  # -- Validate output --------------------------------------------------------
  if (!fs::file_exists(profile_path)) {
    rlang::abort(glue::glue(
      "MetaPhlAn4 profile was not created for sample '{sample_id}'. ",
      "Expected at: {profile_path}"
    ))
  }

  tibble::tibble(
    sample_id        = sample_id,
    profile_path     = as.character(profile_path),
    bowtie2_out_path = as.character(bowtie2_out)
  )
}


#' Parse MetaPhlAn4 profile to tidy format
#'
#' Reads a MetaPhlAn4 output file, filters to species-level clades (those
#' containing `"s__"` but not `"t__"`), and returns a tidy tibble.
#'
#' @param profile_path Character scalar. Path to MetaPhlAn4 TSV output.
#' @param sample_id Character scalar. Sample identifier to include in the
#'   output.
#'
#' @return A [tibble::tibble] with columns: `sample_id`, `clade_name`,
#'   `taxonomy`, `relative_abundance`.
#'
#' @export
parse_metaphlan_profile <- function(profile_path, sample_id) {

  if (!fs::file_exists(profile_path)) {
    rlang::abort(glue::glue("MetaPhlAn profile not found: {profile_path}"))
  }

  # -- Read the profile -------------------------------------------------------
  # MetaPhlAn4 output is tab-separated with a header block starting with '#'.
  # The main columns are: #clade_name, clade_taxid, relative_abundance, ...
  raw_lines <- readr::read_lines(profile_path)

  # Find the header line (starts with #clade_name or #SampleID)
  header_idx <- which(
    stringr::str_detect(raw_lines, "^#clade_name|^#SampleID")
  )

  if (length(header_idx) == 0) {
    # Fallback: skip all comment lines
    comment_lines <- which(stringr::str_detect(raw_lines, "^#"))
    skip_n <- if (length(comment_lines) > 0) max(comment_lines) else 0L
    header_idx <- skip_n
  } else {
    header_idx <- max(header_idx)
  }

  # Read from the header line onward
  profile <- readr::read_tsv(
    I(raw_lines[header_idx:length(raw_lines)]),
    col_types = readr::cols(.default = readr::col_character()),
    comment   = "",
    show_col_types = FALSE
  )

  # Normalise column names (MetaPhlAn uses #clade_name)
  colnames(profile) <- stringr::str_replace(colnames(profile), "^#", "")
  colnames(profile) <- tolower(colnames(profile))

  # Expect at least clade_name and relative_abundance
  if (!"clade_name" %in% colnames(profile)) {
    rlang::abort(glue::glue(
      "Column 'clade_name' not found in MetaPhlAn profile: {profile_path}"
    ))
  }

  # Identify the abundance column (could be relative_abundance or the sample name)
  abund_col <- if ("relative_abundance" %in% colnames(profile)) {
    "relative_abundance"
  } else {
    # MetaPhlAn sometimes uses the sample ID as the column name
    non_meta_cols <- setdiff(colnames(profile), c("clade_name", "clade_taxid",
                                                   "additional_species"))
    if (length(non_meta_cols) >= 1) non_meta_cols[1] else NA_character_
  }

  if (is.na(abund_col)) {
    rlang::abort(glue::glue(
      "Could not identify abundance column in MetaPhlAn profile: {profile_path}"
    ))
  }

  # -- Filter to species level ------------------------------------------------
  species <- profile |>
    dplyr::filter(
      stringr::str_detect(.data$clade_name, "s__"),
      !stringr::str_detect(.data$clade_name, "t__")
    ) |>
    dplyr::mutate(
      relative_abundance = as.numeric(.data[[abund_col]])
    )

  # Extract clean species name: everything after the last "s__"
  species <- species |>
    dplyr::mutate(
      species_name = stringr::str_extract(.data$clade_name, "s__[^|]+$"),
      species_name = stringr::str_replace(.data$species_name, "^s__", "")
    )

  tibble::tibble(
    sample_id          = sample_id,
    clade_name         = species$species_name,
    taxonomy           = species$clade_name,
    relative_abundance = species$relative_abundance
  )
}


#' Merge MetaPhlAn4 profiles across samples
#'
#' Takes the combined output of multiple [parse_metaphlan_profile()] calls and
#' produces both a wide abundance matrix and a long-format tibble.
#'
#' @param parsed_profiles A [tibble::tibble] produced by binding rows from
#'   [parse_metaphlan_profile()]. Must contain columns `sample_id`,
#'   `clade_name`, and `relative_abundance`.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{abundance_matrix}{Wide-format tibble. Rows = species, columns =
#'       samples. Missing species are filled with 0.}
#'     \item{long_format}{Long-format tibble with `sample_id`, `clade_name`,
#'       and `relative_abundance`.}
#'   }
#'
#' @export
merge_metaphlan_profiles <- function(parsed_profiles) {

  required_cols <- c("sample_id", "clade_name", "relative_abundance")
  missing_cols  <- setdiff(required_cols, colnames(parsed_profiles))
  if (length(missing_cols) > 0) {
    rlang::abort(glue::glue(
      "parsed_profiles is missing required columns: ",
      "{paste(missing_cols, collapse = ', ')}"
    ))
  }

  if (nrow(parsed_profiles) == 0) {
    rlang::warn("parsed_profiles has zero rows; returning empty structures.")
    return(list(
      abundance_matrix = tibble::tibble(clade_name = character()),
      long_format      = tibble::tibble(
        sample_id          = character(),
        clade_name         = character(),
        relative_abundance = numeric()
      )
    ))
  }

  # -- Long format: keep only needed columns, ensure numeric ------------------
  long_fmt <- parsed_profiles |>
    dplyr::select(
      "sample_id", "clade_name", "relative_abundance"
    ) |>
    dplyr::mutate(
      relative_abundance = as.numeric(.data$relative_abundance)
    )

  # -- Wide format (abundance matrix) -----------------------------------------
  wide_fmt <- long_fmt |>
    tidyr::pivot_wider(
      names_from  = "sample_id",
      values_from = "relative_abundance",
      values_fill = 0
    )

  list(
    abundance_matrix = wide_fmt,
    long_format      = long_fmt
  )
}
