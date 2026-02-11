# ============================================================================
# step06_instrain.R - inStrain profile & compare wrappers
# ============================================================================

#' Run inStrain profile on a single sample BAM
#'
#' Executes \code{inStrain profile} inside the configured conda environment to
#' generate per-genome microdiversity statistics from a sorted BAM file aligned
#' to a reference genome database.
#'
#' @param sample_id Character scalar. Sample identifier used to name the output
#'   directory.
#' @param bam_path Character scalar. Path to a coordinate-sorted, indexed BAM
#'   file.
#' @param output_dir Character scalar. Base directory under which the inStrain
#'   output folder will be created (e.g. \code{output/06_instrain}).
#' @param config Pipeline config list as returned by
#'   \code{\link{load_config}}.
#' @return A single-row \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{sample_id}{The sample identifier.}
#'     \item{is_profile_dir}{Absolute path to the inStrain output directory.}
#'     \item{genome_info_path}{Absolute path to the \code{genome_info.tsv}
#'       file inside the output directory.}
#'   }
#' @export
run_instrain_profile <- function(sample_id, bam_path, output_dir, config) {

  # -- validate inputs --------------------------------------------------------
  if (!fs::file_exists(bam_path)) {
    rlang::abort(
      glue::glue("BAM file not found for sample '{sample_id}': {bam_path}"),
      class = "strainflow_file_missing"
    )
  }

  # -- resolve config parameters ----------------------------------------------
  ip <- config$instrain_profile
  conda_env      <- config$tools$instrain_env %||% "instrain"
  reference_fasta <- config$reference$fasta
  stb_file       <- config$reference$stb_file
  processes      <- ip$processes   %||% config$execution$workers %||% 4L
  min_read_ani   <- ip$min_read_ani %||% 0.95

  min_mapq       <- ip$min_mapq     %||% 2L
  min_cov        <- ip$min_coverage %||% 5
  min_freq       <- ip$min_freq     %||% 0.05

  # -- build output path ------------------------------------------------------
  is_profile_dir <- fs::path(output_dir, glue::glue("{sample_id}_IS"))
  fs::dir_create(output_dir, recurse = TRUE)

  # -- construct inStrain command string --------------------------------------
  cmd_str <- glue::glue(
    "inStrain profile {bam_path} {reference_fasta}",
    " -o {is_profile_dir}",
    " -p {processes}",
    " -l {min_read_ani}",
    " --min_mapq {min_mapq}",
    " -c {min_cov}",
    " -f {min_freq}",
    " -s {stb_file}",
    " --skip_plot_generation"
  )

  cli::cli_alert_info(
    "Running inStrain profile for sample {.val {sample_id}}"
  )

  run_conda_cmd(
    cmd       = cmd_str,
    conda_env = conda_env,
    echo      = TRUE,
    timeout_sec = config$execution$timeout_sec %||% 14400
  )

  # -- locate genome_info.tsv -------------------------------------------------
  genome_info_path <- fs::path(
    is_profile_dir, "output", "genome_info.tsv"
  )

  if (!fs::file_exists(genome_info_path)) {
    rlang::warn(
      glue::glue(
        "genome_info.tsv not found at expected path for sample ",
        "'{sample_id}': {genome_info_path}"
      )
    )
  }

  tibble::tibble(
    sample_id        = sample_id,
    is_profile_dir   = as.character(fs::path_abs(is_profile_dir)),
    genome_info_path = as.character(fs::path_abs(genome_info_path))
  )
}


#' Parse inStrain genome_info.tsv
#'
#' Reads the \code{genome_info.tsv} file produced by \code{inStrain profile}
#' and returns a tidy tibble with snake_case column names and the originating
#' sample identifier.
#'
#' @param genome_info_path Character scalar. Path to \code{genome_info.tsv}
#'   inside an inStrain output directory.
#' @param sample_id Character scalar. Sample identifier to add as a column.
#' @return A \code{\link[tibble]{tibble}} with columns including:
#'   \code{sample_id}, \code{genome}, \code{coverage}, \code{breadth},
#'   \code{nucl_diversity}, \code{length}, \code{snv_count}, and all other
#'   columns present in the inStrain output (renamed to snake_case).
#' @export
parse_instrain_genome_info <- function(genome_info_path, sample_id) {

  if (!fs::file_exists(genome_info_path)) {
    rlang::abort(
      glue::glue(
        "genome_info.tsv not found: {genome_info_path}"
      ),
      class = "strainflow_file_missing"
    )
  }

  raw <- readr::read_tsv(
    genome_info_path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )

  # -- clean column names to snake_case ---------------------------------------
  clean_names <- names(raw) |>
    tolower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace("^_|_$", "")

  names(raw) <- clean_names

  # -- coerce numeric columns -------------------------------------------------
  numeric_cols <- c(
    "coverage", "breadth", "nucl_diversity", "length",
    "true_scaffolds", "detected_scaffolds", "coverage_median",
    "coverage_std", "coverage_median", "breadth_minc",
    "breadth_expected", "snv_count", "snv_count_unfiltered",
    "divergent_site_count", "snv_per_bp",
    "mean_microdiversity", "median_microdiversity",
    "mean_clonality", "median_clonality",
    "linked_snv_count", "snv_distance_mean",
    "r2_mean", "d_prime_mean",
    "consensus_divergent_sites", "population_divergent_sites",
    "iraw_unfiltered_pairs", "iraw_within_sample_pairs",
    "linked_snvs"
  )

  present_numeric <- intersect(numeric_cols, names(raw))

  result <- raw |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(present_numeric),
        ~ suppressWarnings(as.numeric(.x))
      )
    ) |>
    tibble::as_tibble()

  # -- prepend sample_id ------------------------------------------------------
  result <- dplyr::mutate(result, sample_id = sample_id, .before = 1L)

  result
}


#' Filter genomes by quality thresholds
#'
#' Removes genomes that do not meet minimum breadth and coverage thresholds,
#' ensuring downstream analyses are restricted to well-covered species.
#'
#' @param genome_info A \code{\link[tibble]{tibble}} as returned by
#'   \code{\link{parse_instrain_genome_info}}.
#' @param min_breadth Numeric scalar (0--1). Minimum breadth of coverage.
#'   Default 0.5.
#' @param min_coverage Numeric scalar. Minimum average read depth. Default 5.
#' @return Filtered \code{\link[tibble]{tibble}} retaining only rows that meet
#'   both thresholds.
#' @export
filter_genomes_by_quality <- function(genome_info,
                                      min_breadth = 0.5,
                                      min_coverage = 5) {

  if (!all(c("breadth", "coverage") %in% names(genome_info))) {
    rlang::abort(
      paste0(
        "genome_info must contain 'breadth' and 'coverage' columns. ",
        "Found: ", paste(names(genome_info), collapse = ", ")
      ),
      class = "strainflow_column_missing"
    )
  }

  filtered <- genome_info |>
    dplyr::filter(
      !is.na(.data$breadth),
      !is.na(.data$coverage),
      .data$breadth  >= min_breadth,
      .data$coverage >= min_coverage
    )

  n_removed <- nrow(genome_info) - nrow(filtered)

  cli::cli_alert_info(
    "Genome quality filter: kept {nrow(filtered)}/{nrow(genome_info)} ",
    "genomes (removed {n_removed} below breadth >= {min_breadth} or ",
    "coverage >= {min_coverage})"
  )

  filtered
}


#' Run inStrain compare on a group of profiles
#'
#' Executes \code{inStrain compare} on two or more inStrain profile directories
#' to compute pairwise strain-level comparisons (popANI, conANI, shared SNPs).
#'
#' @param comparison_group_id Character scalar. Identifier for this comparison
#'   group (e.g. a pair_id from the comparison manifest).
#' @param is_profile_dirs Character vector. Paths to inStrain profile output
#'   directories to compare.
#' @param output_dir Character scalar. Directory where the compare output will
#'   be written.
#' @param config Pipeline config list.
#' @return A single-row \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{comparison_group_id}{The group identifier.}
#'     \item{compare_output_dir}{Absolute path to the inStrain compare output
#'       directory.}
#'     \item{comparisons_table_path}{Absolute path to
#'       \code{comparisonsTable.tsv} inside the output.}
#'   }
#' @export
run_instrain_compare <- function(comparison_group_id,
                                 is_profile_dirs,
                                 output_dir,
                                 config) {

  # -- validate inputs --------------------------------------------------------
  if (length(is_profile_dirs) < 2L) {
    rlang::abort(
      glue::glue(
        "inStrain compare requires at least 2 profile directories. ",
        "Got {length(is_profile_dirs)} for group '{comparison_group_id}'."
      ),
      class = "strainflow_invalid_input"
    )
  }

  missing_dirs <- is_profile_dirs[!fs::dir_exists(is_profile_dirs)]
  if (length(missing_dirs) > 0L) {
    rlang::abort(
      c(
        glue::glue(
          "inStrain profile directories not found for group ",
          "'{comparison_group_id}':"
        ),
        missing_dirs
      ),
      class = "strainflow_file_missing"
    )
  }

  # -- resolve config parameters ----------------------------------------------
  ic <- config$instrain_compare
  conda_env  <- config$tools$instrain_env %||% "instrain"
  stb_file   <- config$reference$stb_file
  processes  <- ic$processes  %||% config$execution$workers %||% 4L
  min_cov    <- ic$min_coverage %||% 5
  min_freq   <- ic$min_freq    %||% 0.05

  # -- build output path ------------------------------------------------------
  compare_output_dir <- fs::path(
    output_dir, glue::glue("{comparison_group_id}_compare")
  )
  fs::dir_create(output_dir, recurse = TRUE)

  # -- construct inStrain compare command -------------------------------------
  # Multiple profile dirs are passed as space-separated arguments to -i
  profile_dirs_str <- paste(is_profile_dirs, collapse = " ")

  cmd_str <- glue::glue(
    "inStrain compare",
    " -i {profile_dirs_str}",
    " -o {compare_output_dir}",
    " -p {processes}",
    " -s {stb_file}",
    " -c {min_cov}",
    " -f {min_freq}"
  )

  cli::cli_alert_info(
    "Running inStrain compare for group {.val {comparison_group_id}} ",
    "({length(is_profile_dirs)} profiles)"
  )

  run_conda_cmd(
    cmd       = cmd_str,
    conda_env = conda_env,
    echo      = TRUE,
    timeout_sec = config$execution$timeout_sec %||% 14400
  )

  # -- locate comparisonsTable.tsv --------------------------------------------
  comparisons_table_path <- fs::path(
    compare_output_dir, "output", "comparisonsTable.tsv"
  )

  if (!fs::file_exists(comparisons_table_path)) {
    rlang::warn(
      glue::glue(
        "comparisonsTable.tsv not found at expected path for group ",
        "'{comparison_group_id}': {comparisons_table_path}"
      )
    )
  }

  tibble::tibble(
    comparison_group_id  = comparison_group_id,
    compare_output_dir   = as.character(fs::path_abs(compare_output_dir)),
    comparisons_table_path = as.character(
      fs::path_abs(comparisons_table_path)
    )
  )
}


#' Parse inStrain comparisonsTable.tsv
#'
#' Reads the \code{comparisonsTable.tsv} produced by \code{inStrain compare}
#' and returns a tidy tibble with cleaned column names and appropriate types.
#'
#' @param comparisons_table_path Character scalar. Path to
#'   \code{comparisonsTable.tsv}.
#' @return A \code{\link[tibble]{tibble}} with columns: \code{name1},
#'   \code{name2}, \code{genome}, \code{coverage_overlap},
#'   \code{compared_bases_count}, \code{percent_genome_compared},
#'   \code{consensus_snps}, \code{population_snps}, \code{popani},
#'   \code{conani}, and any additional columns from the inStrain output.
#' @export
parse_instrain_comparisons <- function(comparisons_table_path) {

  if (!fs::file_exists(comparisons_table_path)) {
    rlang::abort(
      glue::glue(
        "comparisonsTable.tsv not found: {comparisons_table_path}"
      ),
      class = "strainflow_file_missing"
    )
  }

  raw <- readr::read_tsv(
    comparisons_table_path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )

  # -- clean column names to snake_case ---------------------------------------
  clean_names <- names(raw) |>
    tolower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace("^_|_$", "")

  names(raw) <- clean_names

  # -- coerce numeric columns -------------------------------------------------
  numeric_cols <- c(
    "coverage_overlap", "compared_bases_count",
    "percent_genome_compared", "consensus_snps", "population_snps",
    "popani", "conani",
    "percent_compared", "length",
    "compared_bases_fraction",
    "mean_coverage_1", "mean_coverage_2"
  )

  present_numeric <- intersect(numeric_cols, names(raw))

  result <- raw |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(present_numeric),
        ~ suppressWarnings(as.numeric(.x))
      )
    ) |>
    tibble::as_tibble()

  result
}
