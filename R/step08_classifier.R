# ============================================================================
# step08_classifier.R - Multi-Signal Transmission Classifier
#
# Combines popANI (SNV signal from inStrain) and APSS (synteny signal from
# SynTracker) into a unified transmission score. This is the core analytical
# innovation: species-specific weighting of orthogonal strain identity signals.
# ============================================================================

# --------------------------------------------------------------------------
# classify_evolutionary_mode
# --------------------------------------------------------------------------

#' Classify species evolutionary mode based on popANI vs APSS variance
#'
#' For each genome (species), computes the variance of SNV-based (popANI)
#' and synteny-based (APSS) strain identity signals across within-individual
#' comparisons (longitudinal positive controls). Species are classified into
#' evolutionary modes that determine how much weight each signal receives
#' in the integrated transmission score.
#'
#' @param instrain_comparisons Tibble from `parse_instrain_comparisons` with
#'   columns: `genome`, `sample_id_1`, `sample_id_2`, `popANI`,
#'   `compared_bases_count`.
#' @param syntracker_results Tibble from `run_syntracker_analysis` with
#'   columns: `genome`, `sample_id_1`, `sample_id_2`, `apss`,
#'   `n_regions_compared`.
#' @param within_individual Tibble of within-individual longitudinal
#'   comparisons (same person, different timepoints) used as positive
#'   controls. Must contain columns `sample_id_1`, `sample_id_2`. If
#'   `NULL`, all comparisons are used (less accurate).
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{genome}{Species/genome identifier.}
#'     \item{snv_variance}{Variance of `(1 - popANI)` across comparisons.}
#'     \item{synteny_variance}{Variance of `(1 - APSS)` across comparisons.}
#'     \item{evolutionary_mode}{One of `"hypermutator"`, `"hyper_recombinator"`,
#'       or `"balanced"`.}
#'     \item{popani_weight}{Weight assigned to popANI in integrated score.}
#'     \item{apss_weight}{Weight assigned to APSS in integrated score.}
#'     \item{n_comparisons}{Number of comparisons used for classification.}
#'   }
#'
#' @export
classify_evolutionary_mode <- function(instrain_comparisons,
                                       syntracker_results,
                                       within_individual = NULL) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(instrain_comparisons, c("genome", "sample_id_1", "sample_id_2", "popANI"),
                 "instrain_comparisons")
  .validate_cols(syntracker_results, c("genome", "sample_id_1", "sample_id_2", "apss"),
                 "syntracker_results")

  # -- Join inStrain and SynTracker by genome + sample pair -------------------
  # Normalize pair ordering to ensure consistent joins
  instrain_norm <- instrain_comparisons |>
    dplyr::mutate(
      pair_key = purrr::map2_chr(
        .data$sample_id_1, .data$sample_id_2,
        ~ paste0(sort(c(.x, .y)), collapse = "__")
      )
    )

  syntracker_norm <- syntracker_results |>
    dplyr::mutate(
      pair_key = purrr::map2_chr(
        .data$sample_id_1, .data$sample_id_2,
        ~ paste0(sort(c(.x, .y)), collapse = "__")
      )
    )

  combined <- instrain_norm |>
    dplyr::select("genome", "pair_key", "popANI") |>
    dplyr::inner_join(
      syntracker_norm |> dplyr::select("genome", "pair_key", "apss"),
      by = c("genome", "pair_key")
    )

  if (nrow(combined) == 0L) {
    rlang::warn(
      c("No overlapping genome + sample pairs between inStrain and SynTracker.",
        "i" = "Returning empty classification table."),
      class = "strainflow_no_overlap"
    )
    return(tibble::tibble(
      genome = character(),
      snv_variance = double(),
      synteny_variance = double(),
      evolutionary_mode = character(),
      popani_weight = double(),
      apss_weight = double(),
      n_comparisons = integer()
    ))
  }

  # -- Filter to within-individual comparisons if provided --------------------
  if (!is.null(within_individual)) {
    .validate_cols(within_individual, c("sample_id_1", "sample_id_2"),
                   "within_individual")
    within_keys <- within_individual |>
      dplyr::mutate(
        pair_key = purrr::map2_chr(
          .data$sample_id_1, .data$sample_id_2,
          ~ paste0(sort(c(.x, .y)), collapse = "__")
        )
      ) |>
      dplyr::pull(.data$pair_key) |>
      unique()

    combined <- combined |>
      dplyr::filter(.data$pair_key %in% .env$within_keys)

    if (nrow(combined) == 0L) {
      rlang::warn(
        c("No within-individual comparisons found after filtering.",
          "i" = "Falling back to all comparisons for evolutionary mode classification."),
        class = "strainflow_no_within_individual"
      )
      # Re-join without filter
      combined <- instrain_norm |>
        dplyr::select("genome", "pair_key", "popANI") |>
        dplyr::inner_join(
          syntracker_norm |> dplyr::select("genome", "pair_key", "apss"),
          by = c("genome", "pair_key")
        )
    }
  }

  # -- Compute per-genome variance of divergence signals ----------------------
  genome_stats <- combined |>
    dplyr::group_by(.data$genome) |>
    dplyr::summarise(
      snv_variance = ifelse(
        dplyr::n() >= 2L,
        stats::var(1 - .data$popANI, na.rm = TRUE),
        NA_real_
      ),
      synteny_variance = ifelse(
        dplyr::n() >= 2L,
        stats::var(1 - .data$apss, na.rm = TRUE),
        NA_real_
      ),
      n_comparisons = dplyr::n(),
      .groups = "drop"
    )

  # -- Classify evolutionary mode ---------------------------------------------
  # Use medians across all species as thresholds
  med_snv <- stats::median(genome_stats$snv_variance, na.rm = TRUE)
  med_syn <- stats::median(genome_stats$synteny_variance, na.rm = TRUE)

  # Guard against degenerate cases where all variances are zero or NA
  if (is.na(med_snv)) med_snv <- 0
  if (is.na(med_syn)) med_syn <- 0

  result <- genome_stats |>
    dplyr::mutate(
      evolutionary_mode = dplyr::case_when(
        is.na(.data$snv_variance) | is.na(.data$synteny_variance) ~ "balanced",
        .data$snv_variance > .env$med_snv & .data$synteny_variance <= .env$med_syn ~ "hypermutator",
        .data$synteny_variance > .env$med_syn & .data$snv_variance <= .env$med_snv ~ "hyper_recombinator",
        TRUE ~ "balanced"
      ),
      popani_weight = dplyr::case_when(
        .data$evolutionary_mode == "hypermutator"       ~ 0.8,
        .data$evolutionary_mode == "hyper_recombinator" ~ 0.2,
        TRUE                                            ~ 0.5
      ),
      apss_weight = dplyr::case_when(
        .data$evolutionary_mode == "hypermutator"       ~ 0.2,
        .data$evolutionary_mode == "hyper_recombinator" ~ 0.8,
        TRUE                                            ~ 0.5
      )
    )

  cli::cli_inform(c(
    "v" = "Evolutionary mode classified for {nrow(result)} genomes.",
    "i" = glue::glue(
      "hypermutator: {sum(result$evolutionary_mode == 'hypermutator')}, ",
      "hyper_recombinator: {sum(result$evolutionary_mode == 'hyper_recombinator')}, ",
      "balanced: {sum(result$evolutionary_mode == 'balanced')}"
    )
  ))

  result
}

# --------------------------------------------------------------------------
# train_species_weights
# --------------------------------------------------------------------------

#' Train species-specific weights using longitudinal positive controls
#'
#' Uses within-individual longitudinal comparisons (same subject, different
#' timepoints) as positive controls and unrelated pairs as negative controls
#' to optimize the weighting of popANI vs APSS for each species. Maximizes
#' AUC for distinguishing same-strain from different-strain comparisons.
#'
#' @param instrain_comparisons Full comparison table with popANI. Must
#'   include columns: `genome`, `sample_id_1`, `sample_id_2`, `popANI`.
#' @param syntracker_results Full APSS table. Must include columns:
#'   `genome`, `sample_id_1`, `sample_id_2`, `apss`.
#' @param sample_metadata Sample metadata with columns: `sample_id`,
#'   `subject_id`, `timepoint_days`.
#' @param evolutionary_modes Output from [classify_evolutionary_mode()].
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{genome}{Species/genome identifier.}
#'     \item{trained_popani_weight}{Optimized weight for popANI.}
#'     \item{trained_apss_weight}{Optimized weight for APSS.}
#'     \item{n_training_pairs}{Total number of training pairs (positive +
#'       negative).}
#'     \item{auc}{Area under ROC curve achieved by the optimal weights.
#'       `NA` if insufficient data for training.}
#'   }
#'
#' @export
train_species_weights <- function(instrain_comparisons,
                                  syntracker_results,
                                  sample_metadata,
                                  evolutionary_modes) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(instrain_comparisons, c("genome", "sample_id_1", "sample_id_2", "popANI"),
                 "instrain_comparisons")
  .validate_cols(syntracker_results, c("genome", "sample_id_1", "sample_id_2", "apss"),
                 "syntracker_results")
  .validate_cols(sample_metadata, c("sample_id", "subject_id", "timepoint_days"),
                 "sample_metadata")
  .validate_cols(evolutionary_modes, c("genome", "popani_weight", "apss_weight"),
                 "evolutionary_modes")

  # -- Annotate comparisons with subject IDs ----------------------------------
  meta_lookup <- sample_metadata |>
    dplyr::select("sample_id", "subject_id") |>
    dplyr::distinct()

  # Normalize pair ordering and join both datasets
  instrain_ann <- instrain_comparisons |>
    dplyr::mutate(
      pair_key = purrr::map2_chr(
        .data$sample_id_1, .data$sample_id_2,
        ~ paste0(sort(c(.x, .y)), collapse = "__")
      )
    ) |>
    dplyr::left_join(
      meta_lookup |> dplyr::rename(subject_1 = "subject_id"),
      by = c("sample_id_1" = "sample_id")
    ) |>
    dplyr::left_join(
      meta_lookup |> dplyr::rename(subject_2 = "subject_id"),
      by = c("sample_id_2" = "sample_id")
    )

  syntracker_ann <- syntracker_results |>
    dplyr::mutate(
      pair_key = purrr::map2_chr(
        .data$sample_id_1, .data$sample_id_2,
        ~ paste0(sort(c(.x, .y)), collapse = "__")
      )
    )

  # Join popANI and APSS
  training_data <- instrain_ann |>
    dplyr::select("genome", "pair_key", "popANI", "subject_1", "subject_2") |>
    dplyr::inner_join(
      syntracker_ann |> dplyr::select("genome", "pair_key", "apss"),
      by = c("genome", "pair_key")
    ) |>
    dplyr::filter(!is.na(.data$popANI), !is.na(.data$apss))

  if (nrow(training_data) == 0L) {
    rlang::warn(
      c("No overlapping data for weight training.",
        "i" = "Falling back to evolutionary mode defaults for all genomes."),
      class = "strainflow_no_training_data"
    )
    return(
      evolutionary_modes |>
        dplyr::select("genome", "popani_weight", "apss_weight") |>
        dplyr::rename(
          trained_popani_weight = "popani_weight",
          trained_apss_weight = "apss_weight"
        ) |>
        dplyr::mutate(
          n_training_pairs = 0L,
          auc = NA_real_
        )
    )
  }

  # -- Label positive (within-individual) and negative (unrelated) controls ---
  training_data <- training_data |>
    dplyr::mutate(
      is_same_strain = dplyr::case_when(
        is.na(.data$subject_1) | is.na(.data$subject_2) ~ NA,
        .data$subject_1 == .data$subject_2 ~ TRUE,
        TRUE ~ FALSE
      )
    ) |>
    dplyr::filter(!is.na(.data$is_same_strain))

  # -- Train per genome -------------------------------------------------------
  genomes_to_train <- training_data |>
    dplyr::group_by(.data$genome) |>
    dplyr::summarise(
      n_positive = sum(.data$is_same_strain),
      n_negative = sum(!.data$is_same_strain),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$n_positive >= 5L, .data$n_negative >= 5L)

  # Weight grid: 0.1 to 0.9 by 0.1
  weight_grid <- seq(0.1, 0.9, by = 0.1)

  trained <- purrr::map_dfr(genomes_to_train$genome, function(g) {
    gdata <- training_data |>
      dplyr::filter(.data$genome == .env$g)

    best_auc <- -1
    best_w <- 0.5

    for (w in weight_grid) {
      gdata_scored <- gdata |>
        dplyr::mutate(score = .env$w * .data$popANI + (1 - .env$w) * .data$apss)

      auc_val <- tryCatch(
        .compute_auc(
          scores = gdata_scored$score,
          labels = gdata_scored$is_same_strain
        ),
        error = function(e) NA_real_
      )

      if (!is.na(auc_val) && auc_val > best_auc) {
        best_auc <- auc_val
        best_w <- w
      }
    }

    n_total <- nrow(gdata)

    tibble::tibble(
      genome = g,
      trained_popani_weight = best_w,
      trained_apss_weight = 1 - best_w,
      n_training_pairs = n_total,
      auc = if (best_auc < 0) NA_real_ else best_auc
    )
  })

  # -- Fall back to evolutionary mode defaults for untrained species ----------
  all_genomes <- unique(c(
    evolutionary_modes$genome,
    instrain_comparisons$genome,
    syntracker_results$genome
  ))

  untrained_genomes <- setdiff(all_genomes, trained$genome)

  fallback <- evolutionary_modes |>
    dplyr::filter(.data$genome %in% .env$untrained_genomes) |>
    dplyr::select("genome", "popani_weight", "apss_weight") |>
    dplyr::rename(
      trained_popani_weight = "popani_weight",
      trained_apss_weight = "apss_weight"
    ) |>
    dplyr::mutate(
      n_training_pairs = 0L,
      auc = NA_real_
    )

  # For genomes not in evolutionary_modes either, use balanced defaults
  still_missing <- setdiff(untrained_genomes, fallback$genome)
  if (length(still_missing) > 0L) {
    extra_fallback <- tibble::tibble(
      genome = still_missing,
      trained_popani_weight = 0.5,
      trained_apss_weight = 0.5,
      n_training_pairs = 0L,
      auc = NA_real_
    )
    fallback <- dplyr::bind_rows(fallback, extra_fallback)
  }

  result <- dplyr::bind_rows(trained, fallback)

  cli::cli_inform(c(
    "v" = "Species weights trained for {nrow(trained)} genomes (AUC-optimized).",
    "i" = glue::glue(
      "{nrow(fallback)} genomes used evolutionary mode defaults (insufficient training data)."
    )
  ))

  result
}

# --------------------------------------------------------------------------
# compute_transmission_score
# --------------------------------------------------------------------------

#' Compute integrated transmission score
#'
#' Combines popANI (SNV-based) and APSS (synteny-based) strain identity
#' signals into a single transmission score using species-specific weights.
#' The score is on a 0--1 scale where higher values indicate more similar
#' strains (higher probability of shared strain identity).
#'
#' @param instrain_comparisons Tibble with popANI values. Must include
#'   columns: `genome`, `sample_id_1`, `sample_id_2`, `popANI`.
#' @param syntracker_results Tibble with APSS values. Must include
#'   columns: `genome`, `sample_id_1`, `sample_id_2`, `apss`.
#' @param species_weights Tibble from [train_species_weights()] or
#'   [classify_evolutionary_mode()]. Must include `genome` and weight
#'   columns (either `trained_popani_weight`/`trained_apss_weight` or
#'   `popani_weight`/`apss_weight`).
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{genome}{Species/genome identifier.}
#'     \item{sample_id_1}{First sample in comparison.}
#'     \item{sample_id_2}{Second sample in comparison.}
#'     \item{popANI}{Population ANI from inStrain (may be `NA`).}
#'     \item{apss}{Average Pairwise Synteny Score from SynTracker (may be
#'       `NA`).}
#'     \item{popani_weight}{Weight applied to popANI.}
#'     \item{apss_weight}{Weight applied to APSS.}
#'     \item{transmission_score}{Integrated score on 0--1 scale.}
#'   }
#'
#' @export
compute_transmission_score <- function(instrain_comparisons,
                                       syntracker_results,
                                       species_weights) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(instrain_comparisons, c("genome", "sample_id_1", "sample_id_2", "popANI"),
                 "instrain_comparisons")
  .validate_cols(syntracker_results, c("genome", "sample_id_1", "sample_id_2", "apss"),
                 "syntracker_results")

  # Detect weight column naming convention
  weight_cols <- colnames(species_weights)
  if ("trained_popani_weight" %in% weight_cols) {
    sw <- species_weights |>
      dplyr::select("genome",
                     popani_weight = "trained_popani_weight",
                     apss_weight = "trained_apss_weight")
  } else if ("popani_weight" %in% weight_cols) {
    sw <- species_weights |>
      dplyr::select("genome", "popani_weight", "apss_weight")
  } else {
    rlang::abort(
      c("species_weights must contain weight columns.",
        "i" = "Expected 'trained_popani_weight'/'trained_apss_weight' or 'popani_weight'/'apss_weight'."),
      class = "strainflow_missing_weights"
    )
  }

  # -- Normalize pair ordering for consistent joins ---------------------------
  instrain_norm <- instrain_comparisons |>
    dplyr::mutate(
      .s1 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[1]),
      .s2 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[2])
    ) |>
    dplyr::mutate(sample_id_1 = .data$.s1, sample_id_2 = .data$.s2) |>
    dplyr::select(-".s1", -".s2")

  syntracker_norm <- syntracker_results |>
    dplyr::mutate(
      .s1 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[1]),
      .s2 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[2])
    ) |>
    dplyr::mutate(sample_id_1 = .data$.s1, sample_id_2 = .data$.s2) |>
    dplyr::select(-".s1", -".s2")

  # -- Full outer join: keep comparisons even if only one signal is available --
  combined <- dplyr::full_join(
    instrain_norm |> dplyr::select("genome", "sample_id_1", "sample_id_2", "popANI"),
    syntracker_norm |> dplyr::select("genome", "sample_id_1", "sample_id_2", "apss"),
    by = c("genome", "sample_id_1", "sample_id_2")
  )

  # -- Join with species weights ----------------------------------------------
  combined <- combined |>
    dplyr::left_join(sw, by = "genome")

  # Apply balanced defaults for genomes without trained weights
  combined <- combined |>
    dplyr::mutate(
      popani_weight = dplyr::coalesce(.data$popani_weight, 0.5),
      apss_weight = dplyr::coalesce(.data$apss_weight, 0.5)
    )

  # -- Compute transmission score ---------------------------------------------
  # Handle missing values: use whichever signal is available
  combined <- combined |>
    dplyr::mutate(
      transmission_score = dplyr::case_when(
        # Both signals available: weighted combination
        !is.na(.data$popANI) & !is.na(.data$apss) ~
          .data$popani_weight * .data$popANI + .data$apss_weight * .data$apss,
        # Only popANI available
        !is.na(.data$popANI) & is.na(.data$apss) ~
          .data$popANI,
        # Only APSS available
        is.na(.data$popANI) & !is.na(.data$apss) ~
          .data$apss,
        # Neither available
        TRUE ~ NA_real_
      )
    )

  result <- combined |>
    dplyr::select(
      "genome", "sample_id_1", "sample_id_2",
      "popANI", "apss",
      "popani_weight", "apss_weight",
      "transmission_score"
    )

  n_both <- sum(!is.na(result$popANI) & !is.na(result$apss))
  n_popani_only <- sum(!is.na(result$popANI) & is.na(result$apss))
  n_apss_only <- sum(is.na(result$popANI) & !is.na(result$apss))

  cli::cli_inform(c(
    "v" = "Transmission scores computed for {nrow(result)} comparisons.",
    "i" = glue::glue(
      "Both signals: {n_both}, popANI only: {n_popani_only}, APSS only: {n_apss_only}"
    )
  ))

  result
}

# --------------------------------------------------------------------------
# call_shared_strains
# --------------------------------------------------------------------------

#' Call shared strains based on transmission score threshold
#'
#' Applies a threshold to integrated transmission scores to call shared
#' (same) strains between sample pairs. Assigns confidence tiers based
#' on data completeness and signal strength.
#'
#' @param transmission_scores Tibble from [compute_transmission_score()].
#' @param threshold Transmission score threshold for calling shared strains
#'   (default `0.99999`). Comparisons with scores at or above this value
#'   are called as shared.
#' @param min_compared Minimum number of compared bases for inStrain data
#'   to be considered reliable (default `50000`).
#' @param min_regions Minimum number of compared regions for SynTracker
#'   data to be considered reliable (default `20`).
#'
#' @return The input tibble with additional columns:
#'   \describe{
#'     \item{is_shared}{Logical. `TRUE` if transmission score >= threshold.}
#'     \item{confidence}{Character. One of `"high"`, `"medium"`, or `"low"`.}
#'   }
#'
#' @export
call_shared_strains <- function(transmission_scores,
                                threshold = 0.99999,
                                min_compared = 50000,
                                min_regions = 20) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(transmission_scores,
                 c("genome", "sample_id_1", "sample_id_2", "transmission_score"),
                 "transmission_scores")

  if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0 || threshold > 1) {
    rlang::abort(
      "threshold must be a single numeric value between 0 and 1.",
      class = "strainflow_invalid_threshold"
    )
  }

  # -- Call shared strains ----------------------------------------------------
  result <- transmission_scores |>
    dplyr::mutate(
      is_shared = !is.na(.data$transmission_score) & .data$transmission_score >= .env$threshold
    )

  # -- Assign confidence tiers ------------------------------------------------
  # Thresholds for individual signals (used for confidence, not for calling)
  popani_threshold <- 0.99999
  apss_threshold <- 0.99

  # Margin around threshold for "marginal" classification
  margin <- (1 - threshold) * 5  # e.g., for 0.99999 threshold, margin ~ 0.00005

  has_popani <- !is.na(result$popANI)
  has_apss <- !is.na(result$apss)
  has_both <- has_popani & has_apss

  popani_strong <- has_popani & result$popANI >= popani_threshold
  apss_strong <- has_apss & result$apss >= apss_threshold

  score_marginal <- !is.na(result$transmission_score) &
    abs(result$transmission_score - threshold) <= margin

  result <- result |>
    dplyr::mutate(
      confidence = dplyr::case_when(
        # No data at all
        is.na(.data$transmission_score) ~ NA_character_,

        # HIGH: both signals available and both individually strong, not marginal
        .env$has_both & .env$popani_strong & .env$apss_strong & !.env$score_marginal ~ "high",

        # MEDIUM: only one signal but strong, or both signals but one is marginal
        (!.env$has_both & (.env$popani_strong | .env$apss_strong)) |
          (.env$has_both & (.env$popani_strong | .env$apss_strong) & .env$score_marginal) ~ "medium",

        # LOW: marginal scores near threshold, or weak signals
        TRUE ~ "low"
      )
    )

  n_shared <- sum(result$is_shared, na.rm = TRUE)
  n_total <- sum(!is.na(result$transmission_score))

  cli::cli_inform(c(
    "v" = "Shared strains called: {n_shared} / {n_total} comparisons at threshold {threshold}.",
    "i" = glue::glue(
      "Confidence: high={sum(result$confidence == 'high', na.rm = TRUE)}, ",
      "medium={sum(result$confidence == 'medium', na.rm = TRUE)}, ",
      "low={sum(result$confidence == 'low', na.rm = TRUE)}"
    )
  ))

  result
}

# ============================================================================
# Internal helper functions
# ============================================================================

#' Validate required columns in a data frame
#'
#' @param df Data frame to check.
#' @param required Character vector of required column names.
#' @param df_name Name of the data frame for error messages.
#' @return Invisible `TRUE` on success; aborts on failure.
#' @noRd
.validate_cols <- function(df, required, df_name) {
  if (!is.data.frame(df)) {
    rlang::abort(
      glue::glue("`{df_name}` must be a data frame or tibble, not {class(df)[1]}."),
      class = "strainflow_invalid_input"
    )
  }
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0L) {
    rlang::abort(
      c(
        glue::glue("`{df_name}` is missing required columns: {paste(missing, collapse = ', ')}"),
        "i" = glue::glue("Required: {paste(required, collapse = ', ')}")
      ),
      class = "strainflow_missing_cols"
    )
  }
  invisible(TRUE)
}

#' Compute AUC using the trapezoidal rule (no external dependencies)
#'
#' Implements the Mann-Whitney U statistic formulation of AUC: the
#' probability that a randomly chosen positive has a higher score than
#' a randomly chosen negative.
#'
#' @param scores Numeric vector of scores.
#' @param labels Logical vector where `TRUE` = positive, `FALSE` = negative.
#' @return Numeric scalar AUC in [0, 1].
#' @noRd
.compute_auc <- function(scores, labels) {
  # Remove NAs
  valid <- !is.na(scores) & !is.na(labels)
  scores <- scores[valid]
  labels <- labels[valid]

  pos_scores <- scores[labels]
  neg_scores <- scores[!labels]

  n_pos <- length(pos_scores)
  n_neg <- length(neg_scores)

  if (n_pos == 0L || n_neg == 0L) {
    return(NA_real_)
  }

  # Mann-Whitney U statistic approach (exact, O(n*m))
  # For computational efficiency with large datasets, use ranks
  all_scores <- c(pos_scores, neg_scores)
  all_labels <- c(rep(TRUE, n_pos), rep(FALSE, n_neg))

  ranks <- rank(all_scores, ties.method = "average")
  pos_rank_sum <- sum(ranks[seq_len(n_pos)])

  # AUC = (R_pos - n_pos*(n_pos+1)/2) / (n_pos * n_neg)
  auc <- (pos_rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

  # Ensure AUC >= 0.5 (flip if classifier is inverted)
  if (auc < 0.5) auc <- 1 - auc

  auc
}
