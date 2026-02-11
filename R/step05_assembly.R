# ---- step05_assembly.R ----
# Metagenome assembly with MEGAHIT.

#' Run MEGAHIT metagenome assembly
#'
#' Assembles non-host metagenomic reads into contigs using MEGAHIT. Contigs
#' shorter than `min_contig_len` (default 1000 bp) are discarded. Assembly
#' statistics (N50, total length, contig count) are computed from the output
#' FASTA.
#'
#' @param sample_id Character scalar. Sample identifier.
#' @param nonhost_r1 Character scalar. Path to non-host forward reads.
#' @param nonhost_r2 Character scalar. Path to non-host reverse reads.
#' @param output_dir Character scalar. Base output directory. MEGAHIT creates
#'   its own subdirectory; this function wraps that behaviour.
#' @param config Named list. Pipeline configuration. Expected elements under
#'   `config$megahit`:
#'   \describe{
#'     \item{threads}{Integer. Number of threads (default 8).}
#'     \item{min_contig_len}{Integer. Minimum contig length to retain
#'       (default 1000).}
#'     \item{memory}{Numeric. Maximum memory fraction (0-1) or absolute bytes
#'       (optional). Passed as `--memory` if provided.}
#'     \item{k_list}{Character. Comma-separated k-mer sizes (optional). Passed
#'       as `--k-list` if provided (e.g. `"21,41,61,81,101"`).}
#'   }
#'
#' @return A single-row [tibble::tibble] with columns: `sample_id`,
#'   `contigs_path`, `n_contigs`, `total_length`, `n50`.
#'
#' @export
run_megahit <- function(sample_id, nonhost_r1, nonhost_r2,
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
  mh_cfg <- config$megahit %||% list()
  threads        <- mh_cfg$threads        %||% 8L
  min_contig_len <- mh_cfg$min_contig_len %||% 1000L
  memory         <- mh_cfg$memory         %||% NULL
  k_list         <- mh_cfg$k_list         %||% NULL

  # -- Build output paths -----------------------------------------------------
  # MEGAHIT requires its output directory to NOT exist, so we use a sample-

  # specific subdirectory that we remove first if it already exists.
  megahit_outdir <- fs::path(output_dir, glue::glue("{sample_id}_megahit"))

  if (fs::dir_exists(megahit_outdir)) {
    cli::cli_alert_warning(
      "Removing existing MEGAHIT output directory: {megahit_outdir}"
    )
    fs::dir_delete(megahit_outdir)
  }

  # -- Build MEGAHIT command --------------------------------------------------
  args <- c(
    "-1", as.character(nonhost_r1),
    "-2", as.character(nonhost_r2),
    "-o", as.character(megahit_outdir),
    "--min-contig-len", as.character(min_contig_len),
    "-t", as.character(threads)
  )

  if (!is.null(memory)) {
    args <- c(args, "--memory", as.character(memory))
  }
  if (!is.null(k_list) && nzchar(k_list)) {
    args <- c(args, "--k-list", as.character(k_list))
  }

  megahit_log <- fs::path(output_dir, glue::glue("{sample_id}_megahit.log"))
  run_cmd(cmd = "megahit", args = args, log_file = as.character(megahit_log))

  # -- Locate final contigs ---------------------------------------------------
  contigs_path <- fs::path(megahit_outdir, "final.contigs.fa")

  if (!fs::file_exists(contigs_path)) {
    rlang::abort(glue::glue(
      "MEGAHIT contigs file not found for sample '{sample_id}'. ",
      "Expected at: {contigs_path}"
    ))
  }

  # -- Compute assembly statistics --------------------------------------------
  stats <- compute_assembly_stats(as.character(contigs_path))

  tibble::tibble(
    sample_id    = sample_id,
    contigs_path = as.character(contigs_path),
    n_contigs    = stats$n_contigs,
    total_length = stats$total_length,
    n50          = stats$n50
  )
}


#' Compute assembly statistics from contigs FASTA
#'
#' Parses a FASTA file to extract contig lengths, then calculates summary
#' statistics including count, total length, N50, and maximum contig length.
#'
#' @param contigs_path Character scalar. Path to a FASTA file containing
#'   assembled contigs.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{n_contigs}{Integer. Number of contigs.}
#'     \item{total_length}{Integer. Sum of all contig lengths.}
#'     \item{n50}{Integer. N50 value -- the minimum contig length such that
#'       contigs of that length or longer cover at least 50\% of the total
#'       assembly.}
#'     \item{max_length}{Integer. Length of the longest contig.}
#'   }
#'
#' @export
compute_assembly_stats <- function(contigs_path) {

  if (!fs::file_exists(contigs_path)) {
    rlang::abort(glue::glue("Contigs FASTA not found: {contigs_path}"))
  }

  # -- Parse FASTA to extract contig lengths ----------------------------------
  # We read the file line-by-line and accumulate sequence lengths.
  # This avoids loading Biostrings for a simple stat computation.
  lines <- readr::read_lines(contigs_path)

  if (length(lines) == 0) {
    rlang::warn(glue::glue("Contigs FASTA is empty: {contigs_path}"))
    return(list(
      n_contigs    = 0L,
      total_length = 0L,
      n50          = NA_integer_,
      max_length   = NA_integer_
    ))
  }

  # Identify header lines
  is_header <- stringr::str_detect(lines, "^>")
  header_indices <- which(is_header)

  if (length(header_indices) == 0) {
    rlang::warn(glue::glue(
      "No FASTA headers found in: {contigs_path}. File may be malformed."
    ))
    return(list(
      n_contigs    = 0L,
      total_length = 0L,
      n50          = NA_integer_,
      max_length   = NA_integer_
    ))
  }

  # Calculate the length of each contig by summing sequence line lengths
  # between consecutive headers.
  n_contigs <- length(header_indices)

  # Define start/end of sequence blocks
  seq_starts <- header_indices + 1L
  seq_ends   <- c(header_indices[-1] - 1L, length(lines))

  contig_lengths <- vapply(seq_len(n_contigs), function(i) {
    if (seq_starts[i] > seq_ends[i]) {
      # Header with no sequence
      return(0L)
    }
    seq_lines <- lines[seq_starts[i]:seq_ends[i]]
    # Remove any whitespace and sum character counts
    sum(nchar(stringr::str_replace_all(seq_lines, "\\s", "")))
  }, integer(1))

  # -- Compute statistics -----------------------------------------------------
  total_length <- sum(contig_lengths)
  max_length   <- max(contig_lengths)

  # N50 calculation: sort lengths descending, cumulative sum, find the length
  # at which 50% of total assembly is covered.
  sorted_lengths <- sort(contig_lengths, decreasing = TRUE)
  cumulative     <- cumsum(sorted_lengths)
  half_total     <- total_length / 2
  n50_idx        <- which(cumulative >= half_total)[1]
  n50            <- sorted_lengths[n50_idx]

  list(
    n_contigs    = as.integer(n_contigs),
    total_length = as.integer(total_length),
    n50          = as.integer(n50),
    max_length   = as.integer(max_length)
  )
}
