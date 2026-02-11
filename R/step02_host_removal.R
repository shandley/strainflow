# ---- step02_host_removal.R ----
# Remove host (human) reads using bowtie2 + samtools.

#' Remove host reads via bowtie2 alignment to human genome
#'
#' Aligns trimmed reads to a human reference genome index with bowtie2
#' (`--very-sensitive`), then extracts unaligned (non-host) read pairs
#' using samtools. The pipeline:
#'
#' ```
#' bowtie2 ... | samtools view -b -f 12 -F 256 |
#'   samtools sort -n | samtools fastq -1 out_r1 -2 out_r2
#' ```
#'
#' Samtools flags:
#' * `-f 12` : both reads in pair unmapped
#' * `-F 256` : exclude secondary alignments
#'
#' @param sample_id Character scalar. Sample identifier.
#' @param trimmed_r1 Character scalar. Path to trimmed forward reads.
#' @param trimmed_r2 Character scalar. Path to trimmed reverse reads.
#' @param output_dir Character scalar. Output directory (created if needed).
#' @param config Named list. Pipeline configuration. Expected elements under
#'   `config$host_removal`:
#'   \describe{
#'     \item{human_index}{Character. Path prefix to the bowtie2 human genome
#'       index.}
#'     \item{threads}{Integer. Number of bowtie2 threads (default 8).}
#'     \item{samtools_threads}{Integer. Threads for samtools sort (default 4).}
#'   }
#'
#' @return A single-row [tibble::tibble] with columns: `sample_id`,
#'   `nonhost_r1`, `nonhost_r2`, `total_reads`, `host_reads`,
#'   `nonhost_reads`, `pct_host`.
#'
#' @export
run_host_removal <- function(sample_id, trimmed_r1, trimmed_r2,
                             output_dir, config) {

  # -- Validate inputs --------------------------------------------------------
  if (!fs::file_exists(trimmed_r1)) {
    rlang::abort(glue::glue("Trimmed R1 not found: {trimmed_r1}"))
  }
  if (!fs::file_exists(trimmed_r2)) {
    rlang::abort(glue::glue("Trimmed R2 not found: {trimmed_r2}"))
  }
  fs::dir_create(output_dir)

  # -- Resolve config ---------------------------------------------------------
  hr_cfg <- config$host_removal %||% list()
  human_index     <- hr_cfg$human_index
  if (is.null(human_index) || !nzchar(human_index)) {
    rlang::abort(
      "config$host_removal$human_index must be set to the bowtie2 human index prefix."
    )
  }
  threads         <- hr_cfg$threads          %||% 8L
  sam_threads     <- hr_cfg$samtools_threads  %||% 4L

  # -- Build output paths -----------------------------------------------------
  nonhost_r1 <- fs::path(output_dir, glue::glue("{sample_id}_nonhost_R1.fastq.gz"))
  nonhost_r2 <- fs::path(output_dir, glue::glue("{sample_id}_nonhost_R2.fastq.gz"))
  bt2_log    <- fs::path(output_dir, glue::glue("{sample_id}_bowtie2_host.log"))

  # -- Build piped command ----------------------------------------------------
  # We redirect bowtie2 stderr to a log file so we can parse alignment stats.
  pipe_cmd <- glue::glue(
    "bowtie2 --very-sensitive ",
    "-x {human_index} ",
    "-1 {trimmed_r1} ",
    "-2 {trimmed_r2} ",
    "--threads {threads} ",
    "2> {bt2_log} ",
    "| samtools view -b -f 12 -F 256 ",
    "| samtools sort -n -@ {sam_threads} ",
    "| samtools fastq ",
    "-1 {nonhost_r1} ",
    "-2 {nonhost_r2} ",
    "-0 /dev/null ",
    "-s /dev/null"
  )

  # Run via bash -c so the pipe is handled by the shell
  host_log <- fs::path(output_dir, glue::glue("{sample_id}_host_removal.log"))
  run_cmd(cmd = "bash", args = c("-c", pipe_cmd), log_file = as.character(host_log))

  # -- Parse bowtie2 log for alignment statistics -----------------------------
  stats <- parse_bowtie2_log(bt2_log, sample_id)

  host_reads    <- stats$aligned_reads
  total_reads   <- stats$total_reads
  nonhost_reads <- total_reads - host_reads
  pct_host      <- if (total_reads > 0) {
    round(host_reads / total_reads * 100, 2)
  } else {
    NA_real_
  }

  tibble::tibble(
    sample_id     = sample_id,
    nonhost_r1    = as.character(nonhost_r1),
    nonhost_r2    = as.character(nonhost_r2),
    total_reads   = total_reads,
    host_reads    = host_reads,
    nonhost_reads = nonhost_reads,
    pct_host      = pct_host
  )
}


#' Parse bowtie2 alignment log
#'
#' Extracts total reads and aligned reads from bowtie2 stderr output.
#'
#' @param log_path Path to saved bowtie2 stderr log.
#' @param sample_id Sample identifier (for error messages).
#' @return Named list with `total_reads` and `aligned_reads`.
#' @keywords internal
parse_bowtie2_log <- function(log_path, sample_id) {

  if (!fs::file_exists(log_path)) {
    rlang::warn(glue::glue(
      "Bowtie2 log not found for sample '{sample_id}': {log_path}. ",
      "Returning NA alignment stats."
    ))
    return(list(total_reads = NA_real_, aligned_reads = NA_real_))
  }

  lines <- readr::read_lines(log_path)

  # Total reads: e.g. "1234567 reads; of these:"

  total_line <- lines[stringr::str_detect(lines, "reads; of these:")]
  total_reads <- if (length(total_line) == 1) {
    as.numeric(stringr::str_extract(total_line, "^\\s*([0-9]+)\\s", group = 1))
  } else {
    NA_real_
  }

  # Overall alignment rate: e.g. "95.23% overall alignment rate"
  rate_line <- lines[stringr::str_detect(lines, "overall alignment rate")]
  if (length(rate_line) == 1 && !is.na(total_reads)) {
    pct <- as.numeric(stringr::str_extract(rate_line, "([0-9.]+)%", group = 1))
    aligned_reads <- round(total_reads * pct / 100)
  } else {
    # Fallback: sum concordant + discordant alignments
    # "12345 (98.76%) aligned concordantly exactly 1 time"
    conc_once <- extract_bt2_count(lines, "aligned concordantly exactly 1 time")
    conc_multi <- extract_bt2_count(lines, "aligned concordantly >1 times")
    disc_once <- extract_bt2_count(lines, "aligned discordantly 1 time")

    aligned_reads <- if (!is.na(conc_once)) {
      sum(conc_once, conc_multi, disc_once, na.rm = TRUE)
    } else {
      NA_real_
    }
  }

  list(total_reads = total_reads, aligned_reads = aligned_reads)
}


#' Extract a count from bowtie2 log lines
#'
#' Searches for a specific pattern in bowtie2 stderr lines and extracts the
#' leading integer.
#'
#' @param lines Character vector of log lines.
#' @param pattern Character scalar. Pattern to match (partial, not regex).
#' @return Numeric scalar or `NA_real_`.
#' @keywords internal
extract_bt2_count <- function(lines, pattern) {
  matched <- lines[stringr::str_detect(lines, stringr::fixed(pattern))]
  if (length(matched) >= 1) {
    as.numeric(stringr::str_extract(matched[1], "^\\s*([0-9]+)", group = 1))
  } else {
    NA_real_
  }
}
