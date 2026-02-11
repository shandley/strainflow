# ---- step04_read_mapping.R ----
# Map reads to reference genomes with bowtie2 + samtools.

#' Map reads to reference genomes with bowtie2
#'
#' Aligns non-host reads to a species/strain reference genome index using
#' bowtie2 in `--very-sensitive` mode, then sorts and indexes the resulting
#' BAM file with samtools.
#'
#' Pipeline:
#' ```
#' bowtie2 -x index -1 r1 -2 r2 --no-unal --very-sensitive -X max_insert |
#'   samtools sort -@ 4 -o out.bam && samtools index out.bam
#' ```
#'
#' @param sample_id Character scalar. Sample identifier.
#' @param nonhost_r1 Character scalar. Path to non-host forward reads.
#' @param nonhost_r2 Character scalar. Path to non-host reverse reads.
#' @param output_dir Character scalar. Output directory (created if needed).
#' @param config Named list. Pipeline configuration. Expected elements under
#'   `config$read_mapping`:
#'   \describe{
#'     \item{species_index}{Character. Path prefix to the bowtie2 reference
#'       genome index.}
#'     \item{threads}{Integer. Bowtie2 alignment threads (default 8).}
#'     \item{samtools_threads}{Integer. Samtools sort threads (default 4).}
#'     \item{max_insert_size}{Integer. Maximum insert size for paired-end
#'       alignment (default 1000).}
#'   }
#'
#' @return A single-row [tibble::tibble] with columns: `sample_id`,
#'   `bam_path`, `bai_path`, `total_reads`, `mapped_reads`, `pct_mapped`.
#'
#' @export
run_read_mapping <- function(sample_id, nonhost_r1, nonhost_r2,
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
  rm_cfg <- config$read_mapping %||% list()

  species_index <- rm_cfg$species_index
  if (is.null(species_index) || !nzchar(species_index)) {
    rlang::abort(
      "config$read_mapping$species_index must be set to the bowtie2 reference index prefix."
    )
  }
  threads         <- rm_cfg$threads          %||% 8L
  sam_threads     <- rm_cfg$samtools_threads  %||% 4L
  max_insert      <- rm_cfg$max_insert_size   %||% 1000L

  # -- Build output paths -----------------------------------------------------
  bam_path <- fs::path(output_dir, glue::glue("{sample_id}_mapped.bam"))
  bai_path <- fs::path(output_dir, glue::glue("{sample_id}_mapped.bam.bai"))
  bt2_log  <- fs::path(output_dir, glue::glue("{sample_id}_bowtie2_mapping.log"))

  # -- Build piped command ----------------------------------------------------
  # Alignment + sort, then index as a follow-up command joined by &&.
  # Bowtie2 stderr is redirected to a log file for stat parsing.
  pipe_cmd <- glue::glue(
    "bowtie2 ",
    "-x {species_index} ",
    "-1 {nonhost_r1} ",
    "-2 {nonhost_r2} ",
    "--threads {threads} ",
    "--no-unal ",
    "--very-sensitive ",
    "-X {max_insert} ",
    "2> {bt2_log} ",
    "| samtools sort -@ {sam_threads} -o {bam_path} ",
    "&& samtools index {bam_path}"
  )

  mapping_log <- fs::path(output_dir, glue::glue("{sample_id}_read_mapping.log"))
  run_cmd(cmd = "bash", args = c("-c", pipe_cmd), log_file = as.character(mapping_log))

  # -- Parse bowtie2 log for mapping statistics -------------------------------
  stats <- parse_mapping_log(bt2_log, sample_id)

  total_reads  <- stats$total_reads
  mapped_reads <- stats$mapped_reads
  pct_mapped   <- if (!is.na(total_reads) && total_reads > 0) {
    round(mapped_reads / total_reads * 100, 2)
  } else {
    NA_real_
  }

  # -- Validate outputs -------------------------------------------------------
  if (!fs::file_exists(bam_path)) {
    rlang::abort(glue::glue(
      "BAM file was not created for sample '{sample_id}': {bam_path}"
    ))
  }
  if (!fs::file_exists(bai_path)) {
    rlang::warn(glue::glue(
      "BAM index was not created for sample '{sample_id}': {bai_path}. ",
      "Some downstream tools may fail."
    ))
  }

  tibble::tibble(
    sample_id    = sample_id,
    bam_path     = as.character(bam_path),
    bai_path     = as.character(bai_path),
    total_reads  = total_reads,
    mapped_reads = mapped_reads,
    pct_mapped   = pct_mapped
  )
}


#' Parse bowtie2 mapping log for alignment statistics
#'
#' Reads the stderr log produced by bowtie2 during reference mapping and
#' extracts total read pairs and the number of aligned reads.
#'
#' @param log_path Character scalar. Path to the bowtie2 stderr log file.
#' @param sample_id Character scalar. Sample identifier (for diagnostics).
#'
#' @return Named list with `total_reads` and `mapped_reads`.
#' @keywords internal
parse_mapping_log <- function(log_path, sample_id) {

  if (!fs::file_exists(log_path)) {
    rlang::warn(glue::glue(
      "Bowtie2 mapping log not found for sample '{sample_id}': {log_path}. ",
      "Returning NA mapping stats."
    ))
    return(list(total_reads = NA_real_, mapped_reads = NA_real_))
  }

  lines <- readr::read_lines(log_path)

  # Total reads: "1234567 reads; of these:"
  total_line <- lines[stringr::str_detect(lines, "reads; of these:")]
  total_reads <- if (length(total_line) >= 1) {
    as.numeric(stringr::str_extract(total_line[1], "^\\s*([0-9]+)\\s", group = 1))
  } else {
    NA_real_
  }

  # Overall alignment rate: "85.23% overall alignment rate"
  rate_line <- lines[stringr::str_detect(lines, "overall alignment rate")]
  if (length(rate_line) >= 1 && !is.na(total_reads)) {
    pct <- as.numeric(
      stringr::str_extract(rate_line[1], "([0-9.]+)%", group = 1)
    )
    mapped_reads <- round(total_reads * pct / 100)
  } else {
    # Fallback: count aligned concordantly + discordantly
    conc_once  <- extract_mapping_count(lines, "aligned concordantly exactly 1 time")
    conc_multi <- extract_mapping_count(lines, "aligned concordantly >1 times")
    disc_once  <- extract_mapping_count(lines, "aligned discordantly 1 time")
    mapped_reads <- if (!is.na(conc_once)) {
      sum(conc_once, conc_multi, disc_once, na.rm = TRUE)
    } else {
      NA_real_
    }
  }

  list(total_reads = total_reads, mapped_reads = mapped_reads)
}


#' Extract a numeric count from bowtie2 log lines
#'
#' Helper to find a line matching a fixed pattern and return the leading integer.
#'
#' @param lines Character vector of log lines.
#' @param pattern Character scalar. Fixed substring to search for.
#' @return Numeric scalar, or `NA_real_` if not found.
#' @keywords internal
extract_mapping_count <- function(lines, pattern) {
  matched <- lines[stringr::str_detect(lines, stringr::fixed(pattern))]
  if (length(matched) >= 1) {
    as.numeric(stringr::str_extract(matched[1], "^\\s*([0-9]+)", group = 1))
  } else {
    NA_real_
  }
}
