# ============================================================================
# step10_longitudinal.R - Longitudinal Strain Dynamics
#
# Tracks the fate of bacterial strains over time in mother-child pairs:
# acquisition timing, persistence/loss, post-transmission divergence,
# and strain replacement events. Uses survival analysis to model strain
# persistence and detects evolutionary divergence from the transmitted
# founding lineage.
# ============================================================================

# --------------------------------------------------------------------------
# build_strain_state_matrix
# --------------------------------------------------------------------------

#' Build strain state matrix tracking presence/absence over time
#'
#' Constructs a longitudinal state matrix recording whether each strain
#' (species) is detected at each timepoint for both mothers and children.
#' Uses within-individual transmission scores to determine whether a
#' strain observed at one timepoint is the "same" strain observed at
#' another.
#'
#' @param shared_strains Tibble from [call_shared_strains()] with all
#'   pairwise comparisons (mother-child, within-mother, within-child).
#'   Must include: `genome`, `sample_id_1`, `sample_id_2`,
#'   `transmission_score`, `is_shared`.
#' @param sample_metadata Sample metadata with: `sample_id`, `subject_id`,
#'   `pair_id`, `role`, `timepoint_days`.
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{pair_id}{Mother-child pair identifier.}
#'     \item{genome}{Species/genome identifier.}
#'     \item{subject_id}{Subject identifier.}
#'     \item{role}{`"mother"` or `"child"`.}
#'     \item{timepoint_days}{Sampling timepoint in days.}
#'     \item{strain_detected}{Logical indicating strain presence.}
#'     \item{transmission_score}{Best transmission score to any other
#'       sample of the same subject for this genome.}
#'   }
#'
#' @export
build_strain_state_matrix <- function(shared_strains, sample_metadata) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(shared_strains,
                 c("genome", "sample_id_1", "sample_id_2",
                   "transmission_score", "is_shared"),
                 "shared_strains")
  .validate_cols(sample_metadata,
                 c("sample_id", "subject_id", "pair_id", "role", "timepoint_days"),
                 "sample_metadata")

  # -- Identify all unique sample-genome combinations -------------------------
  # A genome is "present" in a sample if it appears in any comparison
  # involving that sample (inStrain/SynTracker could map it)
  samples_with_genomes <- dplyr::bind_rows(
    shared_strains |>
      dplyr::select("genome", sample_id = "sample_id_1") |>
      dplyr::distinct(),
    shared_strains |>
      dplyr::select("genome", sample_id = "sample_id_2") |>
      dplyr::distinct()
  ) |>
    dplyr::distinct()

  # -- Annotate with metadata -------------------------------------------------
  meta_lookup <- sample_metadata |>
    dplyr::select("sample_id", "subject_id", "pair_id", "role", "timepoint_days")

  state_matrix <- samples_with_genomes |>
    dplyr::left_join(meta_lookup, by = "sample_id") |>
    dplyr::filter(!is.na(.data$subject_id))

  if (nrow(state_matrix) == 0L) {
    rlang::warn(
      "No samples could be matched to metadata. Returning empty state matrix.",
      class = "strainflow_no_state_data"
    )
    return(tibble::tibble(
      pair_id = character(), genome = character(),
      subject_id = character(), role = character(),
      timepoint_days = double(), strain_detected = logical(),
      transmission_score = double()
    ))
  }

  # -- Determine strain detection per subject x genome x timepoint ------------
  # For within-individual comparisons, strain is "detected" if the
  # genome was profiled in that sample. We use the best transmission
  # score to any other sample of the same subject as a quality indicator.

  # Build a lookup of within-individual transmission scores
  within_scores <- shared_strains |>
    dplyr::left_join(
      meta_lookup |> dplyr::select("sample_id", subject_1 = "subject_id"),
      by = c("sample_id_1" = "sample_id")
    ) |>
    dplyr::left_join(
      meta_lookup |> dplyr::select("sample_id", subject_2 = "subject_id"),
      by = c("sample_id_2" = "sample_id")
    ) |>
    dplyr::filter(
      !is.na(.data$subject_1), !is.na(.data$subject_2),
      .data$subject_1 == .data$subject_2
    )

  # For each sample x genome, get the best within-individual score
  best_within_score_s1 <- within_scores |>
    dplyr::group_by(.data$genome, sample_id = .data$sample_id_1) |>
    dplyr::summarise(
      best_within_score = max(.data$transmission_score, na.rm = TRUE),
      .groups = "drop"
    )

  best_within_score_s2 <- within_scores |>
    dplyr::group_by(.data$genome, sample_id = .data$sample_id_2) |>
    dplyr::summarise(
      best_within_score = max(.data$transmission_score, na.rm = TRUE),
      .groups = "drop"
    )

  best_within_scores <- dplyr::bind_rows(
    best_within_score_s1, best_within_score_s2
  ) |>
    dplyr::group_by(.data$genome, .data$sample_id) |>
    dplyr::summarise(
      best_within_score = max(.data$best_within_score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      # Handle -Inf from empty groups
      best_within_score = dplyr::if_else(
        is.infinite(.data$best_within_score), NA_real_, .data$best_within_score
      )
    )

  # -- Merge scores into state matrix -----------------------------------------
  state_matrix <- state_matrix |>
    dplyr::left_join(
      best_within_scores,
      by = c("genome", "sample_id")
    ) |>
    dplyr::mutate(
      # Strain is detected if it was profiled (genome appeared in comparisons)
      # This is a presence/absence call based on inStrain/SynTracker coverage
      strain_detected = TRUE,
      transmission_score = .data$best_within_score
    ) |>
    dplyr::select(
      "pair_id", "genome", "subject_id", "role",
      "timepoint_days", "strain_detected", "transmission_score"
    )

  # -- Deduplicate (one row per pair x genome x subject x timepoint) ----------
  state_matrix <- state_matrix |>
    dplyr::group_by(
      .data$pair_id, .data$genome, .data$subject_id,
      .data$role, .data$timepoint_days
    ) |>
    dplyr::summarise(
      strain_detected = any(.data$strain_detected),
      transmission_score = max(.data$transmission_score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      transmission_score = dplyr::if_else(
        is.infinite(.data$transmission_score), NA_real_, .data$transmission_score
      )
    )

  cli::cli_inform(c(
    "v" = glue::glue(
      "Strain state matrix built: {nrow(state_matrix)} observations ",
      "({length(unique(state_matrix$genome))} genomes, ",
      "{length(unique(state_matrix$subject_id))} subjects)."
    )
  ))

  state_matrix
}

# --------------------------------------------------------------------------
# detect_acquisition_events
# --------------------------------------------------------------------------

#' Detect acquisition events in children
#'
#' For each species in each child, identifies the first timepoint at
#' which the strain was detected and classifies the acquisition timing
#' relative to birth. Also checks whether the mother carried the same
#' strain at or before the time of acquisition.
#'
#' @param state_matrix Tibble from [build_strain_state_matrix()].
#' @param sample_metadata Sample metadata with: `sample_id`, `subject_id`,
#'   `pair_id`, `role`, `timepoint_days`.
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{pair_id}{Mother-child pair identifier.}
#'     \item{genome}{Species/genome identifier.}
#'     \item{acquisition_timepoint_days}{First timepoint where strain was
#'       detected in child (days).}
#'     \item{acquisition_type}{One of `"birth"` (<= 7 days),
#'       `"delayed"` (> 7 days), or `"unknown"` (only one timepoint).}
#'     \item{maternal_strain_present}{Logical indicating whether the mother
#'       had the strain at or before the child's acquisition timepoint.}
#'     \item{n_child_timepoints}{Number of timepoints sampled for this
#'       child.}
#'   }
#'
#' @export
detect_acquisition_events <- function(state_matrix, sample_metadata) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(state_matrix,
                 c("pair_id", "genome", "subject_id", "role",
                   "timepoint_days", "strain_detected"),
                 "state_matrix")
  .validate_cols(sample_metadata,
                 c("sample_id", "subject_id", "pair_id", "role", "timepoint_days"),
                 "sample_metadata")

  # -- Split into child and mother data ---------------------------------------
  child_data <- state_matrix |>
    dplyr::filter(.data$role == "child", .data$strain_detected)

  mother_data <- state_matrix |>
    dplyr::filter(.data$role == "mother", .data$strain_detected)

  if (nrow(child_data) == 0L) {
    rlang::warn(
      "No strain detections in children. Returning empty acquisition table.",
      class = "strainflow_no_child_detections"
    )
    return(tibble::tibble(
      pair_id = character(), genome = character(),
      acquisition_timepoint_days = double(),
      acquisition_type = character(),
      maternal_strain_present = logical(),
      n_child_timepoints = integer()
    ))
  }

  # -- Find first detection per child x genome --------------------------------
  # Count total child timepoints sampled (whether strain detected or not)
  all_child_timepoints <- state_matrix |>
    dplyr::filter(.data$role == "child") |>
    dplyr::group_by(.data$pair_id, .data$genome, .data$subject_id) |>
    dplyr::summarise(
      n_child_timepoints = dplyr::n_distinct(.data$timepoint_days),
      .groups = "drop"
    )

  first_detection <- child_data |>
    dplyr::group_by(.data$pair_id, .data$genome, .data$subject_id) |>
    dplyr::summarise(
      acquisition_timepoint_days = min(.data$timepoint_days, na.rm = TRUE),
      .groups = "drop"
    )

  # -- Classify acquisition type ----------------------------------------------
  first_detection <- first_detection |>
    dplyr::left_join(
      all_child_timepoints,
      by = c("pair_id", "genome", "subject_id")
    ) |>
    dplyr::mutate(
      acquisition_type = dplyr::case_when(
        .data$n_child_timepoints == 1L ~ "unknown",
        .data$acquisition_timepoint_days <= 7 ~ "birth",
        .data$acquisition_timepoint_days > 7 ~ "delayed",
        TRUE ~ "unknown"
      )
    )

  # -- Check maternal strain presence -----------------------------------------
  # Mother had strain at or before child's acquisition timepoint
  maternal_presence <- mother_data |>
    dplyr::select("pair_id", "genome", mother_timepoint = "timepoint_days") |>
    dplyr::distinct()

  acquisition_events <- first_detection |>
    dplyr::left_join(
      maternal_presence,
      by = c("pair_id", "genome"),
      relationship = "many-to-many"
    ) |>
    dplyr::group_by(
      .data$pair_id, .data$genome, .data$subject_id,
      .data$acquisition_timepoint_days, .data$acquisition_type,
      .data$n_child_timepoints
    ) |>
    dplyr::summarise(
      maternal_strain_present = any(
        !is.na(.data$mother_timepoint) &
          .data$mother_timepoint <= .data$acquisition_timepoint_days
      ),
      .groups = "drop"
    )

  result <- acquisition_events |>
    dplyr::select(
      "pair_id", "genome",
      "acquisition_timepoint_days", "acquisition_type",
      "maternal_strain_present", "n_child_timepoints"
    )

  n_birth <- sum(result$acquisition_type == "birth", na.rm = TRUE)
  n_delayed <- sum(result$acquisition_type == "delayed", na.rm = TRUE)
  n_maternal <- sum(result$maternal_strain_present, na.rm = TRUE)

  cli::cli_inform(c(
    "v" = "Acquisition events detected: {nrow(result)} events across {length(unique(result$pair_id))} pairs.",
    "i" = glue::glue(
      "Birth acquisitions: {n_birth}, Delayed: {n_delayed}, ",
      "Maternal strain present: {n_maternal}/{nrow(result)}"
    )
  ))

  result
}

# --------------------------------------------------------------------------
# fit_persistence_model
# --------------------------------------------------------------------------

#' Fit Kaplan-Meier and Cox PH persistence models
#'
#' Models the persistence (duration of carriage) of transmitted strains
#' in children using survival analysis. Fits both a non-parametric
#' Kaplan-Meier estimator and a Cox proportional hazards model with
#' covariates including species identity, delivery mode, and feeding mode.
#'
#' @param state_matrix Tibble from [build_strain_state_matrix()].
#' @param acquisition_events Tibble from [detect_acquisition_events()].
#' @param pair_metadata Pair metadata with covariates: `pair_id`,
#'   and optionally `delivery_mode`, `feeding_mode`.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{km_fit}{A `survfit` object from [survival::survfit()], or
#'       `NULL` if fitting failed.}
#'     \item{cox_fit}{A `coxph` object from [survival::coxph()], or
#'       `NULL` if fitting failed.}
#'     \item{persistence_summary}{A [tibble::tibble] with columns:
#'       `genome`, `n_acquired`, `n_lost`, `n_censored`,
#'       `median_persistence_days`, `q25_persistence_days`,
#'       `q75_persistence_days`.}
#'   }
#'
#' @export
fit_persistence_model <- function(state_matrix,
                                  acquisition_events,
                                  pair_metadata) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(state_matrix,
                 c("pair_id", "genome", "role", "timepoint_days", "strain_detected"),
                 "state_matrix")
  .validate_cols(acquisition_events,
                 c("pair_id", "genome", "acquisition_timepoint_days"),
                 "acquisition_events")
  .validate_cols(pair_metadata, c("pair_id"), "pair_metadata")

  rlang::check_installed("survival", reason = "to fit persistence models")

  # -- Build survival data ----------------------------------------------------
  # For each acquired strain: time = last detected timepoint - acquisition
  # Event = 1 if strain was lost (not detected at a later timepoint)
  # Event = 0 if still detected at last sample (right-censored)

  child_data <- state_matrix |>
    dplyr::filter(.data$role == "child")

  survival_data <- acquisition_events |>
    dplyr::inner_join(
      child_data |>
        dplyr::select("pair_id", "genome", "timepoint_days", "strain_detected"),
      by = c("pair_id", "genome"),
      relationship = "many-to-many"
    ) |>
    # Keep only timepoints at or after acquisition
    dplyr::filter(.data$timepoint_days >= .data$acquisition_timepoint_days)

  if (nrow(survival_data) == 0L) {
    rlang::warn(
      "No longitudinal data after acquisition events. Cannot fit persistence model.",
      class = "strainflow_no_survival_data"
    )
    return(list(
      km_fit = NULL,
      cox_fit = NULL,
      persistence_summary = tibble::tibble(
        genome = character(), n_acquired = integer(),
        n_lost = integer(), n_censored = integer(),
        median_persistence_days = double(),
        q25_persistence_days = double(),
        q75_persistence_days = double()
      )
    ))
  }

  # For each acquisition event, determine:
  # - Last timepoint where strain was still detected
  # - Whether strain was lost at any subsequent timepoint
  event_data <- survival_data |>
    dplyr::group_by(.data$pair_id, .data$genome, .data$acquisition_timepoint_days) |>
    dplyr::summarise(
      last_detected = max(
        .data$timepoint_days[.data$strain_detected],
        na.rm = TRUE
      ),
      last_sampled = max(.data$timepoint_days, na.rm = TRUE),
      any_timepoints_after_acquisition = sum(
        .data$timepoint_days > .data$acquisition_timepoint_days
      ) > 0L,
      # Strain was lost if there exists a timepoint after the last detection
      lost_after_last_detection = any(
        .data$timepoint_days > max(
          .data$timepoint_days[.data$strain_detected],
          na.rm = TRUE
        ) &
          !.data$strain_detected
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      # Handle -Inf from empty max
      last_detected = dplyr::if_else(
        is.infinite(.data$last_detected),
        .data$acquisition_timepoint_days,
        .data$last_detected
      ),
      # Survival time in days from acquisition
      time = .data$last_detected - .data$acquisition_timepoint_days,
      # Event: 1 = strain lost, 0 = censored (still present at last sample)
      event = dplyr::if_else(
        .data$lost_after_last_detection & .data$any_timepoints_after_acquisition,
        1L, 0L
      ),
      # Ensure positive time for survival analysis
      time = pmax(.data$time, 0.5)
    )

  if (nrow(event_data) < 3L) {
    rlang::warn(
      c(
        glue::glue("Only {nrow(event_data)} events available for survival analysis."),
        "i" = "Need at least 3 events. Returning persistence summary only."
      ),
      class = "strainflow_few_events"
    )
  }

  # -- Join with pair metadata for covariates ---------------------------------
  pair_covariates <- pair_metadata |>
    dplyr::select("pair_id", dplyr::any_of(c("delivery_mode", "feeding_mode")))

  event_data <- event_data |>
    dplyr::left_join(pair_covariates, by = "pair_id")

  # -- Compute per-genome persistence summary ---------------------------------
  persistence_summary <- event_data |>
    dplyr::group_by(.data$genome) |>
    dplyr::summarise(
      n_acquired = dplyr::n(),
      n_lost = sum(.data$event == 1L),
      n_censored = sum(.data$event == 0L),
      median_persistence_days = tryCatch(
        stats::median(.data$time, na.rm = TRUE),
        error = function(e) NA_real_
      ),
      q25_persistence_days = tryCatch(
        stats::quantile(.data$time, 0.25, na.rm = TRUE),
        error = function(e) NA_real_
      ),
      q75_persistence_days = tryCatch(
        stats::quantile(.data$time, 0.75, na.rm = TRUE),
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )

  # -- Fit Kaplan-Meier -------------------------------------------------------
  km_fit <- tryCatch(
    {
      surv_obj <- survival::Surv(time = event_data$time, event = event_data$event)
      survival::survfit(surv_obj ~ 1)
    },
    error = function(e) {
      rlang::warn(
        c("Kaplan-Meier fitting failed.", "x" = conditionMessage(e)),
        class = "strainflow_km_error"
      )
      NULL
    }
  )

  # -- Fit Cox PH model -------------------------------------------------------
  cox_fit <- tryCatch(
    {
      # Build formula dynamically based on available covariates and variation
      cox_terms <- character()

      if (length(unique(event_data$genome)) >= 2L) {
        cox_terms <- c(cox_terms, "genome")
      }
      if ("delivery_mode" %in% colnames(event_data) &&
          length(unique(stats::na.omit(event_data$delivery_mode))) >= 2L) {
        cox_terms <- c(cox_terms, "delivery_mode")
      }
      if ("feeding_mode" %in% colnames(event_data) &&
          length(unique(stats::na.omit(event_data$feeding_mode))) >= 2L) {
        cox_terms <- c(cox_terms, "feeding_mode")
      }

      if (length(cox_terms) == 0L) {
        rlang::inform(
          c("!" = "No covariates with sufficient variation for Cox PH model.",
            "i" = "Returning NULL Cox fit."),
          class = "strainflow_no_cox_covariates"
        )
        NULL
      } else if (nrow(event_data) < 10L || sum(event_data$event) < 3L) {
        rlang::inform(
          c("!" = glue::glue(
            "Insufficient events for Cox PH ({sum(event_data$event)} events, ",
            "{nrow(event_data)} observations)."
          )),
          class = "strainflow_cox_insufficient"
        )
        NULL
      } else {
        cox_formula <- stats::as.formula(
          paste0("survival::Surv(time, event) ~ ", paste(cox_terms, collapse = " + "))
        )
        survival::coxph(cox_formula, data = event_data)
      }
    },
    error = function(e) {
      rlang::warn(
        c("Cox PH model fitting failed.", "x" = conditionMessage(e)),
        class = "strainflow_cox_error"
      )
      NULL
    }
  )

  cli::cli_inform(c(
    "v" = glue::glue(
      "Persistence models fitted: {nrow(event_data)} acquisition events, ",
      "{sum(event_data$event)} loss events."
    ),
    "i" = glue::glue(
      "KM fit: {if (!is.null(km_fit)) 'OK' else 'failed'}, ",
      "Cox fit: {if (!is.null(cox_fit)) 'OK' else 'not fitted'}"
    )
  ))

  list(
    km_fit = km_fit,
    cox_fit = cox_fit,
    persistence_summary = persistence_summary
  )
}

# --------------------------------------------------------------------------
# track_divergence
# --------------------------------------------------------------------------

#' Track post-transmission divergence over time
#'
#' For each transmitted strain, tracks how the mother-child popANI and
#' APSS evolve over time as the child's strain diverges from the maternal
#' founding lineage. Computes both absolute values and deltas from the
#' first post-acquisition timepoint.
#'
#' @param instrain_comparisons Full inStrain comparison table. Must include:
#'   `genome`, `sample_id_1`, `sample_id_2`, `popANI`.
#' @param syntracker_results Full SynTracker results. Must include:
#'   `genome`, `sample_id_1`, `sample_id_2`, `apss`.
#' @param acquisition_events Acquisition events from
#'   [detect_acquisition_events()].
#' @param sample_metadata Sample metadata with: `sample_id`, `subject_id`,
#'   `pair_id`, `role`, `timepoint_days`.
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{pair_id}{Mother-child pair identifier.}
#'     \item{genome}{Species/genome identifier.}
#'     \item{days_since_acquisition}{Days since the strain was first
#'       acquired.}
#'     \item{timepoint_days}{Child's age in days.}
#'     \item{popANI}{Population ANI between mother and child at this
#'       timepoint (may be `NA`).}
#'     \item{apss}{Average pairwise synteny score at this timepoint
#'       (may be `NA`).}
#'     \item{delta_popANI}{Change in popANI from first timepoint.}
#'     \item{delta_apss}{Change in APSS from first timepoint.}
#'   }
#'
#' @export
track_divergence <- function(instrain_comparisons,
                             syntracker_results,
                             acquisition_events,
                             sample_metadata) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(instrain_comparisons,
                 c("genome", "sample_id_1", "sample_id_2", "popANI"),
                 "instrain_comparisons")
  .validate_cols(syntracker_results,
                 c("genome", "sample_id_1", "sample_id_2", "apss"),
                 "syntracker_results")
  .validate_cols(acquisition_events,
                 c("pair_id", "genome", "acquisition_timepoint_days"),
                 "acquisition_events")
  .validate_cols(sample_metadata,
                 c("sample_id", "subject_id", "pair_id", "role", "timepoint_days"),
                 "sample_metadata")

  # -- Identify mother and child samples per pair -----------------------------
  meta_lookup <- sample_metadata |>
    dplyr::select("sample_id", "pair_id", "role", "timepoint_days", "subject_id")

  mother_samples <- meta_lookup |>
    dplyr::filter(.data$role == "mother") |>
    dplyr::select(
      mother_sample = "sample_id",
      "pair_id",
      mother_timepoint = "timepoint_days"
    )

  child_samples <- meta_lookup |>
    dplyr::filter(.data$role == "child") |>
    dplyr::select(
      child_sample = "sample_id",
      "pair_id",
      child_timepoint = "timepoint_days"
    )

  # -- Build mother-child comparison pairs ------------------------------------
  mc_pairs <- mother_samples |>
    dplyr::inner_join(child_samples, by = "pair_id", relationship = "many-to-many")

  if (nrow(mc_pairs) == 0L) {
    rlang::warn(
      "No mother-child sample pairs found in metadata.",
      class = "strainflow_no_mc_pairs"
    )
    return(tibble::tibble(
      pair_id = character(), genome = character(),
      days_since_acquisition = double(), timepoint_days = double(),
      popANI = double(), apss = double(),
      delta_popANI = double(), delta_apss = double()
    ))
  }

  # -- Normalize comparison keys for joining ----------------------------------
  instrain_keyed <- instrain_comparisons |>
    dplyr::mutate(
      .s1 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[1]),
      .s2 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[2])
    ) |>
    dplyr::select(genome = "genome", s1 = ".s1", s2 = ".s2", "popANI")

  syntracker_keyed <- syntracker_results |>
    dplyr::mutate(
      .s1 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[1]),
      .s2 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[2])
    ) |>
    dplyr::select(genome = "genome", s1 = ".s1", s2 = ".s2", "apss")

  # Key the mother-child pairs the same way
  mc_keyed <- mc_pairs |>
    dplyr::mutate(
      s1 = purrr::map2_chr(
        .data$mother_sample, .data$child_sample,
        ~ sort(c(.x, .y))[1]
      ),
      s2 = purrr::map2_chr(
        .data$mother_sample, .data$child_sample,
        ~ sort(c(.x, .y))[2]
      )
    )

  # -- Join with inStrain and SynTracker results ------------------------------
  divergence_raw <- mc_keyed |>
    dplyr::left_join(instrain_keyed, by = c("s1", "s2"), relationship = "many-to-many") |>
    dplyr::left_join(syntracker_keyed, by = c("genome", "s1", "s2"))

  # Filter to genomes that have acquisition events
  divergence_raw <- divergence_raw |>
    dplyr::filter(!is.na(.data$genome)) |>
    dplyr::inner_join(
      acquisition_events |>
        dplyr::select("pair_id", "genome", "acquisition_timepoint_days"),
      by = c("pair_id", "genome")
    )

  if (nrow(divergence_raw) == 0L) {
    rlang::warn(
      "No divergence data found for acquired strains.",
      class = "strainflow_no_divergence_data"
    )
    return(tibble::tibble(
      pair_id = character(), genome = character(),
      days_since_acquisition = double(), timepoint_days = double(),
      popANI = double(), apss = double(),
      delta_popANI = double(), delta_apss = double()
    ))
  }

  # -- Compute divergence relative to acquisition ----------------------------
  divergence_raw <- divergence_raw |>
    dplyr::mutate(
      days_since_acquisition = .data$child_timepoint - .data$acquisition_timepoint_days,
      timepoint_days = .data$child_timepoint
    ) |>
    dplyr::filter(.data$days_since_acquisition >= 0)

  # -- Compute deltas from first post-acquisition timepoint -------------------
  # Get baseline values (earliest child timepoint for each pair x genome)
  baselines <- divergence_raw |>
    dplyr::group_by(.data$pair_id, .data$genome) |>
    dplyr::filter(.data$days_since_acquisition == min(.data$days_since_acquisition)) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::select(
      "pair_id", "genome",
      baseline_popANI = "popANI",
      baseline_apss = "apss"
    )

  result <- divergence_raw |>
    dplyr::left_join(baselines, by = c("pair_id", "genome")) |>
    dplyr::mutate(
      delta_popANI = dplyr::case_when(
        !is.na(.data$popANI) & !is.na(.data$baseline_popANI) ~
          .data$popANI - .data$baseline_popANI,
        TRUE ~ NA_real_
      ),
      delta_apss = dplyr::case_when(
        !is.na(.data$apss) & !is.na(.data$baseline_apss) ~
          .data$apss - .data$baseline_apss,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::select(
      "pair_id", "genome",
      "days_since_acquisition", "timepoint_days",
      "popANI", "apss",
      "delta_popANI", "delta_apss"
    ) |>
    dplyr::arrange(.data$pair_id, .data$genome, .data$days_since_acquisition)

  cli::cli_inform(c(
    "v" = glue::glue(
      "Divergence tracked: {nrow(result)} observations across ",
      "{length(unique(result$genome))} genomes and {length(unique(result$pair_id))} pairs."
    )
  ))

  result
}

# --------------------------------------------------------------------------
# detect_strain_replacement
# --------------------------------------------------------------------------

#' Detect strain replacement events
#'
#' Identifies events where a child's strain at one timepoint is replaced
#' by a different strain at a subsequent timepoint. Detection is based on
#' within-child longitudinal popANI: if consecutive timepoints show
#' popANI below a threshold (indicating different strains), a replacement
#' event is recorded. The new strain is then compared to the mother to
#' determine whether it represents re-transmission or novel acquisition.
#'
#' @param state_matrix Strain state matrix from
#'   [build_strain_state_matrix()].
#' @param instrain_comparisons inStrain results for within-child
#'   longitudinal comparisons. Must include: `genome`, `sample_id_1`,
#'   `sample_id_2`, `popANI`.
#' @param sample_metadata Sample metadata with: `sample_id`, `subject_id`,
#'   `pair_id`, `role`, `timepoint_days`.
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{pair_id}{Mother-child pair identifier.}
#'     \item{genome}{Species/genome identifier.}
#'     \item{replacement_timepoint_days}{Timepoint (days) at which the
#'       new strain was detected.}
#'     \item{previous_timepoint_days}{Timepoint (days) of the old strain.}
#'     \item{within_child_popani}{popANI between consecutive child
#'       timepoints (low = different strains).}
#'     \item{old_strain_popani_to_mother}{popANI of old child strain vs
#'       nearest mother sample (may be `NA`).}
#'     \item{new_strain_popani_to_mother}{popANI of new child strain vs
#'       nearest mother sample (may be `NA`).}
#'   }
#'
#' @export
detect_strain_replacement <- function(state_matrix,
                                      instrain_comparisons,
                                      sample_metadata) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(state_matrix,
                 c("pair_id", "genome", "role", "timepoint_days", "strain_detected"),
                 "state_matrix")
  .validate_cols(instrain_comparisons,
                 c("genome", "sample_id_1", "sample_id_2", "popANI"),
                 "instrain_comparisons")
  .validate_cols(sample_metadata,
                 c("sample_id", "subject_id", "pair_id", "role", "timepoint_days"),
                 "sample_metadata")

  # Replacement detection threshold
  replacement_threshold <- 0.99999

  # -- Get child samples with timepoint info ----------------------------------
  child_meta <- sample_metadata |>
    dplyr::filter(.data$role == "child") |>
    dplyr::select("sample_id", "subject_id", "pair_id", "timepoint_days") |>
    dplyr::arrange(.data$pair_id, .data$subject_id, .data$timepoint_days)

  mother_meta <- sample_metadata |>
    dplyr::filter(.data$role == "mother") |>
    dplyr::select(
      mother_sample = "sample_id",
      "pair_id",
      mother_timepoint = "timepoint_days"
    )

  # -- Build within-child consecutive timepoint pairs -------------------------
  child_consecutive <- child_meta |>
    dplyr::group_by(.data$pair_id, .data$subject_id) |>
    dplyr::arrange(.data$timepoint_days) |>
    dplyr::mutate(
      next_sample_id = dplyr::lead(.data$sample_id),
      next_timepoint = dplyr::lead(.data$timepoint_days)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$next_sample_id)) |>
    dplyr::rename(
      prev_sample_id = "sample_id",
      prev_timepoint = "timepoint_days"
    )

  if (nrow(child_consecutive) == 0L) {
    rlang::warn(
      "No consecutive child timepoint pairs found. Cannot detect replacements.",
      class = "strainflow_no_consecutive"
    )
    return(tibble::tibble(
      pair_id = character(), genome = character(),
      replacement_timepoint_days = double(),
      previous_timepoint_days = double(),
      within_child_popani = double(),
      old_strain_popani_to_mother = double(),
      new_strain_popani_to_mother = double()
    ))
  }

  # -- Normalize inStrain comparison keys for lookup --------------------------
  instrain_keyed <- instrain_comparisons |>
    dplyr::mutate(
      s1 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[1]),
      s2 = purrr::map2_chr(.data$sample_id_1, .data$sample_id_2, ~ sort(c(.x, .y))[2])
    ) |>
    dplyr::select("genome", "s1", "s2", "popANI")

  # -- Lookup within-child popANI for consecutive timepoints ------------------
  child_consecutive <- child_consecutive |>
    dplyr::mutate(
      s1 = purrr::map2_chr(
        .data$prev_sample_id, .data$next_sample_id,
        ~ sort(c(.x, .y))[1]
      ),
      s2 = purrr::map2_chr(
        .data$prev_sample_id, .data$next_sample_id,
        ~ sort(c(.x, .y))[2]
      )
    )

  # Cross with genomes that are present in state matrix for these children
  child_genomes <- state_matrix |>
    dplyr::filter(.data$role == "child", .data$strain_detected) |>
    dplyr::select("pair_id", "genome") |>
    dplyr::distinct()

  consecutive_with_genomes <- child_consecutive |>
    dplyr::inner_join(child_genomes, by = "pair_id", relationship = "many-to-many")

  # Join with within-child popANI
  consecutive_scored <- consecutive_with_genomes |>
    dplyr::left_join(instrain_keyed, by = c("genome", "s1", "s2"))

  # -- Detect replacements: popANI below threshold ----------------------------
  replacements <- consecutive_scored |>
    dplyr::filter(
      !is.na(.data$popANI),
      .data$popANI < .env$replacement_threshold
    )

  if (nrow(replacements) == 0L) {
    cli::cli_inform(c(
      "v" = "No strain replacement events detected."
    ))
    return(tibble::tibble(
      pair_id = character(), genome = character(),
      replacement_timepoint_days = double(),
      previous_timepoint_days = double(),
      within_child_popani = double(),
      old_strain_popani_to_mother = double(),
      new_strain_popani_to_mother = double()
    ))
  }

  # -- Get popANI of old and new child strains vs nearest mother sample -------
  # For each replacement, find the closest mother sample to the timepoint
  .get_mother_popani <- function(child_sample_id, pair_id_val, genome_val,
                                 child_tp, mother_meta_df, instrain_keyed_df) {
    # Find mother samples for this pair
    mothers <- mother_meta_df |>
      dplyr::filter(.data$pair_id == .env$pair_id_val)

    if (nrow(mothers) == 0L) return(NA_real_)

    # Pick closest mother sample by timepoint
    mothers <- mothers |>
      dplyr::mutate(
        tp_diff = abs(.data$mother_timepoint - .env$child_tp)
      ) |>
      dplyr::arrange(.data$tp_diff) |>
      dplyr::slice_head(n = 1L)

    mother_sid <- mothers$mother_sample[1]

    # Lookup popANI
    key_s1 <- sort(c(child_sample_id, mother_sid))[1]
    key_s2 <- sort(c(child_sample_id, mother_sid))[2]

    match <- instrain_keyed_df |>
      dplyr::filter(
        .data$genome == .env$genome_val,
        .data$s1 == .env$key_s1,
        .data$s2 == .env$key_s2
      )

    if (nrow(match) == 0L) return(NA_real_)
    match$popANI[1]
  }

  # Compute mother comparisons for old and new strains
  replacements <- replacements |>
    dplyr::mutate(
      old_strain_popani_to_mother = purrr::pmap_dbl(
        list(
          child_sample_id = .data$prev_sample_id,
          pair_id_val = .data$pair_id,
          genome_val = .data$genome,
          child_tp = .data$prev_timepoint
        ),
        ~ .get_mother_popani(
          ..1, ..2, ..3, ..4,
          mother_meta_df = mother_meta,
          instrain_keyed_df = instrain_keyed
        )
      ),
      new_strain_popani_to_mother = purrr::pmap_dbl(
        list(
          child_sample_id = .data$next_sample_id,
          pair_id_val = .data$pair_id,
          genome_val = .data$genome,
          child_tp = .data$next_timepoint
        ),
        ~ .get_mother_popani(
          ..1, ..2, ..3, ..4,
          mother_meta_df = mother_meta,
          instrain_keyed_df = instrain_keyed
        )
      )
    )

  result <- replacements |>
    dplyr::select(
      "pair_id", "genome",
      replacement_timepoint_days = "next_timepoint",
      previous_timepoint_days = "prev_timepoint",
      within_child_popani = "popANI",
      "old_strain_popani_to_mother",
      "new_strain_popani_to_mother"
    ) |>
    dplyr::arrange(.data$pair_id, .data$genome, .data$replacement_timepoint_days)

  # Classify replacement source
  n_from_mother <- sum(
    !is.na(result$new_strain_popani_to_mother) &
      result$new_strain_popani_to_mother >= replacement_threshold,
    na.rm = TRUE
  )
  n_novel <- sum(
    !is.na(result$new_strain_popani_to_mother) &
      result$new_strain_popani_to_mother < replacement_threshold,
    na.rm = TRUE
  )

  cli::cli_inform(c(
    "v" = "Strain replacements detected: {nrow(result)} events.",
    "i" = glue::glue(
      "Re-transmission from mother: {n_from_mother}, ",
      "Novel strain: {n_novel}, ",
      "Unknown source: {nrow(result) - n_from_mother - n_novel}"
    )
  ))

  result
}
