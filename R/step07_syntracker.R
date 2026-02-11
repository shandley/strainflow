# ============================================================================
# step07_syntracker.R - SynTracker-style APSS computation (native R)
#
# Reimplements the core Average Pairwise Synteny Score (APSS) algorithm from
# SynTracker entirely in R using DECIPHER and Biostrings, without calling the
# external SynTracker Python tool.
# ============================================================================

#' Fragment reference genome into query regions for synteny analysis
#'
#' Splits each scaffold of a reference genome FASTA into non-overlapping
#' regions of fixed length separated by a constant spacing. These regions
#' serve as query anchors for detecting synteny rearrangements among
#' metagenomic assemblies.
#'
#' @param reference_fasta Character scalar. Path to the reference genome FASTA
#'   file.
#' @param region_length Integer scalar. Length of each query region in base
#'   pairs. Default 1000.
#' @param region_spacing Integer scalar. Spacing between the start positions
#'   of consecutive regions in base pairs. Default 4000.
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{region_id}{Unique region identifier
#'       (\code{<scaffold>_<start>_<end>}).}
#'     \item{genome}{Name of the FASTA file (basename without extension).}
#'     \item{scaffold}{Scaffold / contig name from the FASTA header.}
#'     \item{start}{1-based start coordinate.}
#'     \item{end}{1-based end coordinate.}
#'     \item{sequence}{Nucleotide sequence of the region as a character
#'       string.}
#'   }
#' @export
fragment_reference <- function(reference_fasta,
                               region_length = 1000L,
                               region_spacing = 4000L) {

  if (!fs::file_exists(reference_fasta)) {
    rlang::abort(
      glue::glue("Reference FASTA not found: {reference_fasta}"),
      class = "strainflow_file_missing"
    )
  }

  rlang::check_installed("Biostrings", reason = "to read FASTA sequences")

  seqs <- Biostrings::readDNAStringSet(reference_fasta)
  genome_name <- tools::file_path_sans_ext(basename(reference_fasta))
  scaffold_names <- names(seqs)
  scaffold_lengths <- Biostrings::width(seqs)

  # -- generate regions for each scaffold ------------------------------------
  regions_list <- vector("list", length(scaffold_names))

  for (i in seq_along(scaffold_names)) {
    scaf_name   <- scaffold_names[i]
    scaf_length <- scaffold_lengths[i]
    scaf_seq    <- seqs[[i]]

    # Start positions: 1, 1 + spacing, 1 + 2*spacing, ...
    starts <- seq(1L, scaf_length, by = as.integer(region_spacing))

    # Keep only regions that fit entirely within the scaffold
    ends <- starts + as.integer(region_length) - 1L
    valid <- ends <= scaf_length
    starts <- starts[valid]
    ends   <- ends[valid]

    if (length(starts) == 0L) next

    # Extract sequences
    region_seqs <- vapply(
      seq_along(starts),
      function(j) {
        as.character(
          Biostrings::subseq(scaf_seq, start = starts[j], end = ends[j])
        )
      },
      character(1)
    )

    region_ids <- glue::glue("{scaf_name}_{starts}_{ends}")

    regions_list[[i]] <- tibble::tibble(
      region_id = as.character(region_ids),
      genome    = genome_name,
      scaffold  = scaf_name,
      start     = starts,
      end       = ends,
      sequence  = region_seqs
    )
  }

  result <- dplyr::bind_rows(regions_list)

  cli::cli_alert_info(
    "Fragmented {.val {genome_name}}: {nrow(result)} regions ",
    "({region_length} bp, {region_spacing} bp spacing) across ",
    "{length(scaffold_names)} scaffold(s)"
  )

  result
}


#' Find homologous regions in a metagenomic assembly via BLAST
#'
#' BLASTs the query regions (from \code{\link{fragment_reference}}) against a
#' metagenomic assembly to identify homologous loci. High-identity hits that
#' cover a sufficient fraction of the query region are retained, and the
#' target sequence plus 2 kb of flanking context is extracted for downstream
#' synteny analysis.
#'
#' @param query_regions A \code{\link[tibble]{tibble}} as returned by
#'   \code{\link{fragment_reference}}.
#' @param assembly_fasta Character scalar. Path to the metagenomic assembly
#'   FASTA file.
#' @param min_identity Numeric scalar (0--1). Minimum BLAST percent identity
#'   for a hit to be retained. Default 0.97.
#' @param min_coverage Numeric scalar (0--1). Minimum fraction of the query
#'   region length covered by the alignment. Default 0.70.
#' @param config Pipeline config list. Used to locate the \code{blastn}
#'   executable and set resource limits. Default empty list.
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{region_id}{Query region identifier.}
#'     \item{genome}{Genome name from the query regions.}
#'     \item{target_scaffold}{Subject scaffold name in the assembly.}
#'     \item{target_start}{1-based start of the BLAST hit on the subject.}
#'     \item{target_end}{1-based end of the BLAST hit on the subject.}
#'     \item{identity}{Percent identity of the alignment.}
#'     \item{query_coverage}{Fraction of query length covered.}
#'     \item{target_sequence_with_flanks}{Subject sequence including 2 kb
#'       upstream and downstream flanking sequence.}
#'   }
#' @export
find_homologs <- function(query_regions,
                          assembly_fasta,
                          min_identity = 0.97,
                          min_coverage = 0.70,
                          config = list()) {

  if (!fs::file_exists(assembly_fasta)) {
    rlang::abort(
      glue::glue("Assembly FASTA not found: {assembly_fasta}"),
      class = "strainflow_file_missing"
    )
  }

  if (nrow(query_regions) == 0L) {
    cli::cli_alert_warning("No query regions provided; returning empty tibble.")
    return(tibble::tibble(
      region_id                 = character(),
      genome                    = character(),
      target_scaffold           = character(),
      target_start              = integer(),
      target_end                = integer(),
      identity                  = double(),
      query_coverage            = double(),
      target_sequence_with_flanks = character()
    ))
  }

  rlang::check_installed("Biostrings", reason = "to extract flanking sequences")

  # -- write query regions to temp FASTA --------------------------------------
  query_fasta <- fs::path(tempdir(), "syntracker_query.fa")
  query_lines <- paste0(">", query_regions$region_id, "\n", query_regions$sequence)
  writeLines(query_lines, query_fasta)

  # -- build BLAST command ----------------------------------------------------
  blast_outfmt <- "6 qseqid sseqid pident length qlen qstart qend sstart send evalue"
  perc_identity_arg <- as.character(round(min_identity * 100, 2))

  blast_result <- tryCatch(
    run_cmd(
      cmd  = "blastn",
      args = c(
        "-query",         as.character(query_fasta),
        "-subject",       assembly_fasta,
        "-outfmt",        blast_outfmt,
        "-perc_identity", perc_identity_arg,
        "-max_target_seqs", "1",
        "-evalue",        "1e-10"
      ),
      echo = FALSE,
      spinner = FALSE,
      timeout_sec = config$execution$timeout_sec %||% 3600
    ),
    error = function(e) {
      rlang::abort(
        c(
          "BLAST search failed.",
          "i" = glue::glue("Query: {query_fasta}"),
          "i" = glue::glue("Subject: {assembly_fasta}"),
          "x" = conditionMessage(e)
        ),
        class = "strainflow_blast_error",
        parent = e
      )
    }
  )

  # -- clean up temp file -----------------------------------------------------
  fs::file_delete(query_fasta)

  # -- parse BLAST output -----------------------------------------------------
  if (nchar(trimws(blast_result$stdout)) == 0L) {
    cli::cli_alert_warning(
      "BLAST returned no hits for assembly {.file {assembly_fasta}}."
    )
    return(tibble::tibble(
      region_id                 = character(),
      genome                    = character(),
      target_scaffold           = character(),
      target_start              = integer(),
      target_end                = integer(),
      identity                  = double(),
      query_coverage            = double(),
      target_sequence_with_flanks = character()
    ))
  }

  hits <- readr::read_tsv(
    I(blast_result$stdout),
    col_names = c(
      "qseqid", "sseqid", "pident", "length",
      "qlen", "qstart", "qend", "sstart", "send", "evalue"
    ),
    col_types = readr::cols(
      qseqid = readr::col_character(),
      sseqid = readr::col_character(),
      pident = readr::col_double(),
      length = readr::col_integer(),
      qlen   = readr::col_integer(),
      qstart = readr::col_integer(),
      qend   = readr::col_integer(),
      sstart = readr::col_integer(),
      send   = readr::col_integer(),
      evalue = readr::col_double()
    ),
    show_col_types = FALSE
  )

  # -- filter by query coverage -----------------------------------------------
  hits <- hits |>
    dplyr::mutate(
      query_coverage = .data$length / .data$qlen,
      identity       = .data$pident / 100
    ) |>
    dplyr::filter(
      .data$query_coverage >= min_coverage,
      .data$identity       >= min_identity
    )

  if (nrow(hits) == 0L) {
    cli::cli_alert_warning(
      "No BLAST hits passed filters for assembly {.file {assembly_fasta}}."
    )
    return(tibble::tibble(
      region_id                 = character(),
      genome                    = character(),
      target_scaffold           = character(),
      target_start              = integer(),
      target_end                = integer(),
      identity                  = double(),
      query_coverage            = double(),
      target_sequence_with_flanks = character()
    ))
  }

  # -- keep best hit per query region -----------------------------------------
  hits <- hits |>
    dplyr::group_by(.data$qseqid) |>
    dplyr::arrange(dplyr::desc(.data$identity), .by_group = TRUE) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup()

  # -- normalise subject coordinates (handle reverse strand) ------------------
  hits <- hits |>
    dplyr::mutate(
      target_start_norm = pmin(.data$sstart, .data$send),
      target_end_norm   = pmax(.data$sstart, .data$send)
    )

  # -- extract target sequences with 2 kb flanking ----------------------------
  assembly_seqs <- Biostrings::readDNAStringSet(assembly_fasta)
  assembly_widths <- stats::setNames(
    Biostrings::width(assembly_seqs),
    names(assembly_seqs)
  )

  flank_bp <- 2000L

  target_seqs <- vapply(
    seq_len(nrow(hits)),
    function(j) {
      scaf <- hits$sseqid[j]

      if (!scaf %in% names(assembly_seqs)) {
        return(NA_character_)
      }

      scaf_len   <- assembly_widths[[scaf]]
      ext_start  <- max(1L, hits$target_start_norm[j] - flank_bp)
      ext_end    <- min(scaf_len, hits$target_end_norm[j] + flank_bp)

      as.character(
        Biostrings::subseq(assembly_seqs[[scaf]], start = ext_start, end = ext_end)
      )
    },
    character(1)
  )

  # -- join genome info from query_regions ------------------------------------
  region_genome_map <- query_regions |>
    dplyr::select("region_id", "genome") |>
    dplyr::distinct()

  result <- tibble::tibble(
    region_id                 = hits$qseqid,
    target_scaffold           = hits$sseqid,
    target_start              = hits$target_start_norm,
    target_end                = hits$target_end_norm,
    identity                  = hits$identity,
    query_coverage            = hits$query_coverage,
    target_sequence_with_flanks = target_seqs
  ) |>
    dplyr::left_join(region_genome_map, by = "region_id") |>
    dplyr::select(
      "region_id", "genome", "target_scaffold", "target_start",
      "target_end", "identity", "query_coverage",
      "target_sequence_with_flanks"
    ) |>
    dplyr::filter(!is.na(.data$target_sequence_with_flanks))

  cli::cli_alert_info(
    "Found {nrow(result)} homologous regions in {.file {assembly_fasta}}"
  )

  result
}


#' Compute synteny score between two homologous regions
#'
#' Aligns two nucleotide sequences using \code{DECIPHER::AlignSeqs()}, then
#' identifies synteny blocks (contiguous stretches of aligned bases) and
#' computes a synteny score following the SynTracker formula:
#'
#' \deqn{synScore = 1 + \log_{10}\!\bigl(\frac{\text{overlap\_length}}
#'   {\text{shorter\_length} \times n\_blocks}\bigr)}{synScore = 1 +
#'   log10(overlap_length / (shorter_length * n_blocks))}
#'
#' A single contiguous alignment block gives a perfect synteny score of 1.0
#' (no rearrangements). Multiple blocks reduce the score, reflecting
#' structural divergence.
#'
#' @param seq1 Character scalar. First nucleotide sequence.
#' @param seq2 Character scalar. Second nucleotide sequence.
#' @return Numeric scalar. Synteny score in the range (0, 1], where 1.0
#'   indicates perfect synteny.
#' @export
compute_synteny_score <- function(seq1, seq2) {

  # -- guard: empty or missing sequences --------------------------------------
  if (is.na(seq1) || is.na(seq2) ||
      nchar(seq1) == 0L || nchar(seq2) == 0L) {
    return(NA_real_)
  }

  rlang::check_installed("Biostrings", reason = "to create DNAStringSet")
  rlang::check_installed("DECIPHER", reason = "to align sequences")

  # -- align the two sequences ------------------------------------------------
  aligned <- tryCatch(
    {
      dna <- Biostrings::DNAStringSet(c(seq1 = seq1, seq2 = seq2))
      DECIPHER::AlignSeqs(dna, verbose = FALSE)
    },
    error = function(e) {
      rlang::warn(
        c(
          "DECIPHER::AlignSeqs failed; returning NA synteny score.",
          "x" = conditionMessage(e)
        )
      )
      return(NULL)
    }
  )

  if (is.null(aligned)) {
    return(NA_real_)
  }

  # -- convert alignment to character matrix ----------------------------------
  aln_mat <- as.matrix(aligned)  # rows = sequences, cols = alignment positions

  seq1_chars <- aln_mat[1, ]
  seq2_chars <- aln_mat[2, ]

  # Positions where both sequences have a non-gap character
  both_present <- (seq1_chars != "-") & (seq2_chars != "-")

  overlap_length <- sum(both_present)

  if (overlap_length == 0L) {
    return(NA_real_)
  }

  # -- identify synteny blocks ------------------------------------------------
  # A synteny block is a maximal contiguous run of alignment columns where
  # both sequences have non-gap characters. A gap in either sequence breaks
  # the block.
  #
  # We detect block boundaries by finding transitions from FALSE to TRUE in

  # the both_present vector.
  n_blocks <- 0L
  in_block <- FALSE

  for (k in seq_along(both_present)) {
    if (both_present[k]) {
      if (!in_block) {
        n_blocks <- n_blocks + 1L
        in_block <- TRUE
      }
    } else {
      in_block <- FALSE
    }
  }

  # -- compute synteny score --------------------------------------------------
  shorter_length <- min(nchar(seq1), nchar(seq2))

  if (shorter_length == 0L || n_blocks == 0L) {
    return(NA_real_)
  }

  # Single block = perfect synteny
  if (n_blocks == 1L) {
    return(1.0)
  }

  # SynTracker formula: synScore = 1 + log10(overlap / (shorter * n_blocks))
  ratio <- overlap_length / (shorter_length * n_blocks)

  # Clamp ratio to (0, 1] to keep score in a meaningful range
  ratio <- min(ratio, 1.0)

  if (ratio <= 0) {
    return(NA_real_)
  }

  syn_score <- 1.0 + log10(ratio)

  # Score should be in (0, 1]; clamp at 0 if log term drives it negative

  syn_score <- max(syn_score, 0.0)

  syn_score
}


#' Compute Average Pairwise Synteny Score (APSS) between two samples
#'
#' For a single species/genome, computes the mean synteny score across shared
#' genomic regions between two samples' metagenomic assemblies. This is the
#' core SynTracker metric for detecting structural variation indicative of
#' strain transmission.
#'
#' @param sample1_homologs A \code{\link[tibble]{tibble}} of homologous
#'   regions for sample 1, as returned by \code{\link{find_homologs}}.
#' @param sample2_homologs A \code{\link[tibble]{tibble}} of homologous
#'   regions for sample 2, as returned by \code{\link{find_homologs}}.
#' @param n_regions Integer scalar. Maximum number of shared regions to
#'   subsample for scoring. If fewer shared regions are available, all are
#'   used. Default 100.
#' @param seed Integer scalar. Random seed for reproducible subsampling.
#'   Default 42.
#' @return Numeric scalar. The APSS value (typically 0--1), or \code{NA} if
#'   no shared regions are available.
#' @export
compute_apss <- function(sample1_homologs,
                         sample2_homologs,
                         n_regions = 100L,
                         seed = 42L) {

  # -- find regions present in both samples -----------------------------------
  shared_regions <- dplyr::inner_join(
    sample1_homologs |>
      dplyr::select(
        "region_id",
        seq1 = "target_sequence_with_flanks"
      ),
    sample2_homologs |>
      dplyr::select(
        "region_id",
        seq2 = "target_sequence_with_flanks"
      ),
    by = "region_id"
  )

  # Remove rows where either sequence is NA
  shared_regions <- shared_regions |>
    dplyr::filter(!is.na(.data$seq1), !is.na(.data$seq2))

  n_shared <- nrow(shared_regions)

  if (n_shared == 0L) {
    cli::cli_alert_warning("No shared regions between the two samples.")
    return(NA_real_)
  }

  # -- subsample if more than n_regions shared regions are available ----------
  if (n_shared > n_regions) {
    set.seed(seed)
    sample_idx <- sample.int(n_shared, size = n_regions)
    shared_regions <- shared_regions[sample_idx, ]
  }

  cli::cli_alert_info(
    "Computing APSS over {nrow(shared_regions)} shared regions ",
    "(of {n_shared} total)."
  )

  # -- compute synteny score for each shared region ---------------------------
  scores <- vapply(
    seq_len(nrow(shared_regions)),
    function(j) {
      compute_synteny_score(
        shared_regions$seq1[j],
        shared_regions$seq2[j]
      )
    },
    numeric(1)
  )

  # -- return mean, ignoring failed alignments --------------------------------
  valid_scores <- scores[!is.na(scores)]

  if (length(valid_scores) == 0L) {
    cli::cli_alert_warning(
      "All synteny score computations returned NA."
    )
    return(NA_real_)
  }

  mean(valid_scores)
}


#' Run full SynTracker analysis for all sample pairs and genomes
#'
#' Orchestrates the complete SynTracker-style APSS analysis:
#' \enumerate{
#'   \item Fragments the reference genome into query regions.
#'   \item Finds homologous regions in each sample's metagenomic assembly via
#'     BLAST.
#'   \item Computes pairwise APSS for every sample pair defined in the
#'     comparison manifest.
#' }
#'
#' @param assembly_results A \code{\link[tibble]{tibble}} with at least
#'   columns \code{sample_id} and \code{contigs_path} (paths to assembled
#'   contigs FASTA files from step 05).
#' @param reference_fasta Character scalar. Path to the reference genome
#'   FASTA.
#' @param sample_metadata A \code{\link[tibble]{tibble}} of sample metadata
#'   (must include \code{sample_id}).
#' @param pair_metadata A \code{\link[tibble]{tibble}} of pairwise
#'   comparisons to perform, with columns \code{pair_id},
#'   \code{sample_id_1}, and \code{sample_id_2}.
#' @param config Pipeline config list.
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{genome}{Reference genome / species name.}
#'     \item{sample_id_1}{First sample in the pair.}
#'     \item{sample_id_2}{Second sample in the pair.}
#'     \item{pair_id}{Pair identifier from the manifest.}
#'     \item{n_regions_compared}{Number of shared regions used for APSS.}
#'     \item{apss}{Average Pairwise Synteny Score.}
#'   }
#' @export
run_syntracker_analysis <- function(assembly_results,
                                    reference_fasta,
                                    sample_metadata,
                                    pair_metadata,
                                    config) {

  # -- validate inputs --------------------------------------------------------
  required_assembly_cols <- c("sample_id", "contigs_path")
  missing_cols <- setdiff(required_assembly_cols, names(assembly_results))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "assembly_results is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  required_pair_cols <- c("pair_id", "sample_id_1", "sample_id_2")
  missing_pair_cols <- setdiff(required_pair_cols, names(pair_metadata))
  if (length(missing_pair_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "pair_metadata is missing required columns: ",
        "{paste(missing_pair_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  # -- resolve syntracker config parameters -----------------------------------
  st_cfg <- config$syntracker
  region_length  <- st_cfg$region_length  %||% 1000L
  region_spacing <- st_cfg$region_spacing %||% 4000L
  min_identity   <- st_cfg$min_identity   %||% 0.97
  min_coverage   <- st_cfg$min_coverage   %||% 0.70
  n_regions      <- st_cfg$n_regions      %||% 100L
  seed           <- st_cfg$seed           %||% 42L

  # -- step 1: fragment reference genome --------------------------------------
  cli::cli_h2("SynTracker: Fragmenting reference genome")
  query_regions <- fragment_reference(
    reference_fasta = reference_fasta,
    region_length   = as.integer(region_length),
    region_spacing  = as.integer(region_spacing)
  )

  # Get unique genome names from the fragments
  genome_names <- unique(query_regions$genome)

  # -- step 2: find homologs in each sample's assembly ------------------------
  cli::cli_h2("SynTracker: Finding homologs in assemblies")

  # Build a named list: sample_id -> homologs tibble
  sample_ids <- assembly_results$sample_id
  homologs_by_sample <- stats::setNames(
    vector("list", length(sample_ids)),
    sample_ids
  )

  for (sid in sample_ids) {
    contigs_path <- assembly_results |>
      dplyr::filter(.data$sample_id == sid) |>
      dplyr::pull("contigs_path")

    if (length(contigs_path) == 0L || !fs::file_exists(contigs_path)) {
      cli::cli_alert_warning(
        "Assembly not found for sample {.val {sid}}; skipping."
      )
      homologs_by_sample[[sid]] <- tibble::tibble(
        region_id                 = character(),
        genome                    = character(),
        target_scaffold           = character(),
        target_start              = integer(),
        target_end                = integer(),
        identity                  = double(),
        query_coverage            = double(),
        target_sequence_with_flanks = character()
      )
      next
    }

    cli::cli_alert_info("Finding homologs for sample {.val {sid}}")

    homologs_by_sample[[sid]] <- tryCatch(
      find_homologs(
        query_regions  = query_regions,
        assembly_fasta = contigs_path,
        min_identity   = min_identity,
        min_coverage   = min_coverage,
        config         = config
      ),
      error = function(e) {
        cli::cli_alert_danger(
          "Homolog search failed for sample {.val {sid}}: {conditionMessage(e)}"
        )
        tibble::tibble(
          region_id                 = character(),
          genome                    = character(),
          target_scaffold           = character(),
          target_start              = integer(),
          target_end                = integer(),
          identity                  = double(),
          query_coverage            = double(),
          target_sequence_with_flanks = character()
        )
      }
    )
  }

  # -- step 3: compute APSS for each sample pair and genome -------------------
  cli::cli_h2("SynTracker: Computing pairwise APSS")

  results_list <- vector("list", nrow(pair_metadata) * length(genome_names))
  result_idx <- 0L

  for (g in genome_names) {
    for (p in seq_len(nrow(pair_metadata))) {
      result_idx <- result_idx + 1L

      pid  <- pair_metadata$pair_id[p]
      sid1 <- pair_metadata$sample_id_1[p]
      sid2 <- pair_metadata$sample_id_2[p]

      # Get homologs for this genome in each sample
      s1_homologs <- homologs_by_sample[[sid1]]
      s2_homologs <- homologs_by_sample[[sid2]]

      if (is.null(s1_homologs) || is.null(s2_homologs)) {
        results_list[[result_idx]] <- tibble::tibble(
          genome             = g,
          sample_id_1        = sid1,
          sample_id_2        = sid2,
          pair_id            = pid,
          n_regions_compared = 0L,
          apss               = NA_real_
        )
        next
      }

      # Filter to current genome
      s1_genome <- s1_homologs |> dplyr::filter(.data$genome == g)
      s2_genome <- s2_homologs |> dplyr::filter(.data$genome == g)

      # Count shared regions for reporting
      shared_region_ids <- intersect(s1_genome$region_id, s2_genome$region_id)
      n_shared <- length(shared_region_ids)

      if (n_shared == 0L) {
        results_list[[result_idx]] <- tibble::tibble(
          genome             = g,
          sample_id_1        = sid1,
          sample_id_2        = sid2,
          pair_id            = pid,
          n_regions_compared = 0L,
          apss               = NA_real_
        )
        next
      }

      cli::cli_alert_info(
        "APSS: genome={.val {g}}, pair={.val {pid}} ",
        "({n_shared} shared regions)"
      )

      apss_value <- tryCatch(
        compute_apss(
          sample1_homologs = s1_genome,
          sample2_homologs = s2_genome,
          n_regions        = as.integer(n_regions),
          seed             = as.integer(seed)
        ),
        error = function(e) {
          cli::cli_alert_danger(
            "APSS computation failed for genome={.val {g}}, pair={.val {pid}}: ",
            "{conditionMessage(e)}"
          )
          NA_real_
        }
      )

      n_compared <- min(n_shared, n_regions)

      results_list[[result_idx]] <- tibble::tibble(
        genome             = g,
        sample_id_1        = sid1,
        sample_id_2        = sid2,
        pair_id            = pid,
        n_regions_compared = as.integer(n_compared),
        apss               = apss_value
      )
    }
  }

  # -- combine results --------------------------------------------------------
  results <- dplyr::bind_rows(results_list)

  cli::cli_alert_success(
    "SynTracker analysis complete: {nrow(results)} genome-pair comparisons"
  )

  results
}
