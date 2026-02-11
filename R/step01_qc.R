# ---- step01_qc.R ----
# Quality control with fastp for paired-end shotgun metagenomics reads.

#' Run fastp quality control on paired-end reads
#'
#' Trims adapters, filters low-quality reads, and generates QC reports using
#' fastp.
#'
#' @param sample_id Character scalar. Sample identifier used to name output
#'   files.
#' @param fastq_r1 Character scalar. Path to forward reads (R1) FASTQ file
#'   (may be gzipped).
#' @param fastq_r2 Character scalar. Path to reverse reads (R2) FASTQ file
#'   (may be gzipped).
#' @param output_dir Character scalar. Directory where trimmed reads and reports
#'   will be written. Created if it does not exist.
#' @param config Named list. Pipeline configuration. Expected elements under
#'   `config$fastp`:
#'   \describe{
#'     \item{threads}{Integer. Number of threads (default 4).}
#'     \item{length_required}{Integer. Minimum read length after trimming
#'       (default 50).}
#'     \item{qualified_quality_phred}{Integer. Phred quality threshold
#'       (default 20).}
#'     \item{cut_window_size}{Integer. Sliding window size for cutting
#'       (default 4).}
#'     \item{cut_mean_quality}{Integer. Mean quality threshold for the sliding
#'       window (default 20).}
#'   }
#'
#' @return A single-row [tibble::tibble] with columns: `sample_id`,
#'   `trimmed_r1`, `trimmed_r2`, `json_report`, `html_report`,
#'   `reads_before`, `reads_after`, `pct_passing`.
#'
#' @export
run_fastp <- function(sample_id, fastq_r1, fastq_r2, output_dir, config) {


  # -- Validate inputs --------------------------------------------------------
  if (!fs::file_exists(fastq_r1)) {
    rlang::abort(glue::glue("Forward reads file not found: {fastq_r1}"))
  }
  if (!fs::file_exists(fastq_r2)) {
    rlang::abort(glue::glue("Reverse reads file not found: {fastq_r2}"))
  }
  fs::dir_create(output_dir)


  # -- Resolve config ---------------------------------------------------------
  fp_cfg <- config$fastp %||% list()
  threads        <- fp_cfg$threads                %||% 4L
  len_req        <- fp_cfg$length_required        %||% 50L
  qual_phred     <- fp_cfg$qualified_quality_phred %||% 20L
  cut_win        <- fp_cfg$cut_window_size        %||% 4L
  cut_mean       <- fp_cfg$cut_mean_quality       %||% 20L

  # -- Build output paths -----------------------------------------------------
  trimmed_r1  <- fs::path(output_dir, glue::glue("{sample_id}_trimmed_R1.fastq.gz"))
  trimmed_r2  <- fs::path(output_dir, glue::glue("{sample_id}_trimmed_R2.fastq.gz"))
  json_report <- fs::path(output_dir, glue::glue("{sample_id}_fastp.json"))
  html_report <- fs::path(output_dir, glue::glue("{sample_id}_fastp.html"))

  # -- Assemble fastp command -------------------------------------------------
  args <- c(
    "--in1",  fastq_r1,
    "--in2",  fastq_r2,
    "--out1", trimmed_r1,
    "--out2", trimmed_r2,
    "--json", json_report,
    "--html", html_report,
    "--thread",                  as.character(threads),
    "--length_required",         as.character(len_req),
    "--qualified_quality_phred", as.character(qual_phred),
    "--cut_front",
    "--cut_tail",
    "--cut_window_size",         as.character(cut_win),
    "--cut_mean_quality",        as.character(cut_mean),
    "--detect_adapter_for_pe"
  )

  fastp_log <- fs::path(output_dir, glue::glue("{sample_id}_fastp.log"))
  run_cmd(cmd = "fastp", args = args, log_file = as.character(fastp_log))

  # -- Parse JSON report ------------------------------------------------------
  if (!fs::file_exists(json_report)) {
    rlang::abort(glue::glue(
      "fastp JSON report was not created for sample '{sample_id}'. ",
      "Expected at: {json_report}"
    ))
  }

  report   <- jsonlite::fromJSON(json_report, simplifyVector = TRUE)
  reads_before <- report$summary$before_filtering$total_reads
  reads_after  <- report$summary$after_filtering$total_reads

  if (is.null(reads_before) || is.null(reads_after)) {
    rlang::warn(glue::glue(
      "Could not parse read counts from fastp JSON for sample '{sample_id}'."
    ))
    reads_before <- NA_real_
    reads_after  <- NA_real_
  }

  pct_passing <- if (!is.na(reads_before) && reads_before > 0) {
    round(reads_after / reads_before * 100, 2)
  } else {
    NA_real_
  }

  # -- Return tibble ----------------------------------------------------------
  tibble::tibble(
    sample_id   = sample_id,
    trimmed_r1  = as.character(trimmed_r1),
    trimmed_r2  = as.character(trimmed_r2),
    json_report = as.character(json_report),
    html_report = as.character(html_report),
    reads_before = as.numeric(reads_before),
    reads_after  = as.numeric(reads_after),
    pct_passing  = pct_passing
  )
}


#' Aggregate QC statistics across samples
#'
#' Takes the combined output of multiple [run_fastp()] calls and computes
#' summary statistics: total reads before/after, overall pass rate, and
#' per-sample quantiles.
#'
#' @param qc_results A [tibble::tibble] produced by binding rows from multiple
#'   [run_fastp()] calls. Must contain columns `sample_id`, `reads_before`,
#'   `reads_after`, and `pct_passing`.
#'
#' @return A [tibble::tibble] with summary statistics:
#'   \describe{
#'     \item{n_samples}{Number of samples.}
#'     \item{total_reads_before}{Sum of reads before trimming.}
#'     \item{total_reads_after}{Sum of reads after trimming.}
#'     \item{overall_pct_passing}{Overall percentage of reads passing QC.}
#'     \item{median_pct_passing}{Median per-sample pass rate.}
#'     \item{min_pct_passing}{Minimum per-sample pass rate.}
#'     \item{max_pct_passing}{Maximum per-sample pass rate.}
#'     \item{mean_reads_before}{Mean reads per sample before trimming.}
#'     \item{mean_reads_after}{Mean reads per sample after trimming.}
#'   }
#'
#' @export
aggregate_qc_stats <- function(qc_results) {

  required_cols <- c("sample_id", "reads_before", "reads_after", "pct_passing")
  missing_cols  <- setdiff(required_cols, colnames(qc_results))
  if (length(missing_cols) > 0) {
    rlang::abort(glue::glue(
      "qc_results is missing required columns: {paste(missing_cols, collapse = ', ')}"
    ))
  }

  if (nrow(qc_results) == 0) {
    rlang::warn("qc_results has zero rows; returning empty summary.")
    return(tibble::tibble(
      n_samples            = 0L,
      total_reads_before   = NA_real_,
      total_reads_after    = NA_real_,
      overall_pct_passing  = NA_real_,
      median_pct_passing   = NA_real_,
      min_pct_passing      = NA_real_,
      max_pct_passing      = NA_real_,
      mean_reads_before    = NA_real_,
      mean_reads_after     = NA_real_
    ))
  }

  total_before <- sum(qc_results$reads_before, na.rm = TRUE)
  total_after  <- sum(qc_results$reads_after,  na.rm = TRUE)
  overall_pct  <- if (total_before > 0) {
    round(total_after / total_before * 100, 2)
  } else {
    NA_real_
  }

  tibble::tibble(
    n_samples            = nrow(qc_results),
    total_reads_before   = total_before,
    total_reads_after    = total_after,
    overall_pct_passing  = overall_pct,
    median_pct_passing   = round(stats::median(qc_results$pct_passing, na.rm = TRUE), 2),
    min_pct_passing      = round(min(qc_results$pct_passing, na.rm = TRUE), 2),
    max_pct_passing      = round(max(qc_results$pct_passing, na.rm = TRUE), 2),
    mean_reads_before    = round(mean(qc_results$reads_before, na.rm = TRUE), 0),
    mean_reads_after     = round(mean(qc_results$reads_after,  na.rm = TRUE), 0)
  )
}
