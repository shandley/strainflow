# ============================================================================
# step09_transmission_model.R - Transmission vs Environmental Acquisition
#
# Bayesian framework for distinguishing true mother-to-child strain
# transmission from coincidental strain sharing due to environmental
# acquisition. Uses background sharing rates from unrelated pairs as
# a null model and incorporates temporal priors based on birth proximity.
# ============================================================================

# --------------------------------------------------------------------------
# estimate_strain_prevalence
# --------------------------------------------------------------------------

#' Estimate background strain prevalence from unrelated control pairs
#'
#' Computes the background rate of strain sharing among unrelated
#' (negative control) pairs for each species. This serves as the null
#' expectation: the probability that two unrelated individuals share the
#' same strain by chance (e.g., due to common environmental sources).
#'
#' @param shared_strains Tibble from [call_shared_strains()] containing all
#'   comparisons including negative controls. Must include columns: `genome`,
#'   `sample_id_1`, `sample_id_2`, `is_shared`.
#' @param sample_metadata Sample metadata with columns: `sample_id`,
#'   `subject_id`, `pair_id`.
#' @param pair_metadata Pair metadata with columns: `pair_id`,
#'   `mother_subject_id`, `child_subject_id`.
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{genome}{Species/genome identifier.}
#'     \item{n_unrelated_pairs_tested}{Number of unrelated pair comparisons
#'       evaluated for this genome.}
#'     \item{n_unrelated_shared}{Number of unrelated pairs sharing this
#'       strain.}
#'     \item{background_sharing_rate}{Proportion of unrelated pairs sharing
#'       the strain.}
#'     \item{se}{Standard error of the sharing rate estimate.}
#'   }
#'
#' @export
estimate_strain_prevalence <- function(shared_strains,
                                       sample_metadata,
                                       pair_metadata) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(shared_strains,
                 c("genome", "sample_id_1", "sample_id_2", "is_shared"),
                 "shared_strains")
  .validate_cols(sample_metadata, c("sample_id", "subject_id", "pair_id"),
                 "sample_metadata")
  .validate_cols(pair_metadata, c("pair_id", "mother_subject_id", "child_subject_id"),
                 "pair_metadata")

  # -- Identify unrelated comparisons -----------------------------------------
  # Build set of related subject pairs from pair_metadata
  related_pairs <- pair_metadata |>
    dplyr::select("mother_subject_id", "child_subject_id") |>
    dplyr::mutate(
      related_key = purrr::map2_chr(
        .data$mother_subject_id, .data$child_subject_id,
        ~ paste0(sort(c(.x, .y)), collapse = "__")
      )
    ) |>
    dplyr::pull(.data$related_key) |>
    unique()

  # Also include within-individual pairs as non-negative-control
  within_individual_keys <- sample_metadata |>
    dplyr::select("sample_id", "subject_id") |>
    dplyr::distinct()

  # Annotate shared_strains with subject IDs
  annotated <- shared_strains |>
    dplyr::left_join(
      within_individual_keys |> dplyr::rename(subject_1 = "subject_id"),
      by = c("sample_id_1" = "sample_id")
    ) |>
    dplyr::left_join(
      within_individual_keys |> dplyr::rename(subject_2 = "subject_id"),
      by = c("sample_id_2" = "sample_id")
    )

  # Classify comparison type
  annotated <- annotated |>
    dplyr::mutate(
      subject_pair_key = purrr::map2_chr(
        .data$subject_1, .data$subject_2,
        ~ if (is.na(.x) || is.na(.y)) NA_character_
          else paste0(sort(c(.x, .y)), collapse = "__")
      ),
      is_within_individual = !is.na(.data$subject_1) & !is.na(.data$subject_2) &
        .data$subject_1 == .data$subject_2,
      is_related = !is.na(.data$subject_pair_key) &
        .data$subject_pair_key %in% .env$related_pairs,
      is_unrelated = !.data$is_within_individual & !.data$is_related &
        !is.na(.data$subject_1) & !is.na(.data$subject_2)
    )

  # Filter to unrelated comparisons only
  unrelated <- annotated |>
    dplyr::filter(.data$is_unrelated)

  if (nrow(unrelated) == 0L) {
    rlang::warn(
      c("No unrelated (negative control) comparisons found.",
        "i" = "Background prevalence cannot be estimated. Returning zeros."),
      class = "strainflow_no_neg_controls"
    )
    genomes <- unique(shared_strains$genome)
    return(tibble::tibble(
      genome = genomes,
      n_unrelated_pairs_tested = 0L,
      n_unrelated_shared = 0L,
      background_sharing_rate = 0,
      se = NA_real_
    ))
  }

  # -- Compute per-genome background sharing rate -----------------------------
  prevalence <- unrelated |>
    dplyr::group_by(.data$genome) |>
    dplyr::summarise(
      n_unrelated_pairs_tested = dplyr::n(),
      n_unrelated_shared = sum(.data$is_shared, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      background_sharing_rate = .data$n_unrelated_shared / .data$n_unrelated_pairs_tested,
      se = sqrt(.data$background_sharing_rate * (1 - .data$background_sharing_rate) /
                  .data$n_unrelated_pairs_tested)
    )

  # Add genomes with no unrelated data
  missing_genomes <- setdiff(unique(shared_strains$genome), prevalence$genome)
  if (length(missing_genomes) > 0L) {
    prevalence <- dplyr::bind_rows(
      prevalence,
      tibble::tibble(
        genome = missing_genomes,
        n_unrelated_pairs_tested = 0L,
        n_unrelated_shared = 0L,
        background_sharing_rate = 0,
        se = NA_real_
      )
    )
  }

  cli::cli_inform(c(
    "v" = "Background prevalence estimated for {nrow(prevalence)} genomes.",
    "i" = glue::glue(
      "Median background sharing rate: ",
      "{round(stats::median(prevalence$background_sharing_rate, na.rm = TRUE), 5)}"
    )
  ))

  prevalence
}

# --------------------------------------------------------------------------
# compute_transmission_odds
# --------------------------------------------------------------------------

#' Compute Bayesian transmission odds
#'
#' For each shared strain event in related (mother-child) pairs, computes
#' the odds that the observed sharing is due to vertical transmission
#' versus coincidental environmental acquisition. Incorporates a temporal
#' prior that decays with the child's age (birth events are more likely
#' to represent transmission than later acquisitions).
#'
#' @param shared_strains Tibble with sharing calls for related pairs. Must
#'   include columns: `genome`, `sample_id_1`, `sample_id_2`,
#'   `transmission_score`, `is_shared`.
#' @param prevalence Tibble from [estimate_strain_prevalence()].
#' @param sample_metadata Sample metadata with columns: `sample_id`,
#'   `subject_id`, `pair_id`, `role`, `timepoint_days`.
#' @param prior_weight Prior weight on transmission at birth (default `0.5`
#'   for uninformative prior). Higher values encode stronger belief that
#'   sharing events represent true transmission.
#'
#' @return A [tibble::tibble] with columns:
#'   \describe{
#'     \item{genome}{Species/genome identifier.}
#'     \item{pair_id}{Mother-child pair identifier.}
#'     \item{sample_id_mother}{Mother sample in comparison.}
#'     \item{sample_id_child}{Child sample in comparison.}
#'     \item{timepoint_child}{Child's age in days at sampling.}
#'     \item{observed_sharing}{Logical, whether strain was called as
#'       shared.}
#'     \item{transmission_score}{Integrated transmission score.}
#'     \item{background_rate}{Background sharing rate for this species.}
#'     \item{prior}{Temporal prior probability of transmission.}
#'     \item{transmission_odds}{Odds ratio favoring transmission over
#'       chance.}
#'     \item{log_odds}{log10 of transmission_odds.}
#'   }
#'
#' @export
compute_transmission_odds <- function(shared_strains,
                                      prevalence,
                                      sample_metadata,
                                      prior_weight = 0.5) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(shared_strains,
                 c("genome", "sample_id_1", "sample_id_2",
                   "transmission_score", "is_shared"),
                 "shared_strains")
  .validate_cols(prevalence, c("genome", "background_sharing_rate"), "prevalence")
  .validate_cols(sample_metadata,
                 c("sample_id", "subject_id", "pair_id", "role", "timepoint_days"),
                 "sample_metadata")

  if (!is.numeric(prior_weight) || length(prior_weight) != 1L ||
      prior_weight <= 0 || prior_weight >= 1) {
    rlang::abort(
      "prior_weight must be a single numeric value in (0, 1).",
      class = "strainflow_invalid_prior"
    )
  }

  # -- Annotate comparisons with mother/child roles ---------------------------
  meta_lookup <- sample_metadata |>
    dplyr::select("sample_id", "subject_id", "pair_id", "role", "timepoint_days")

  annotated <- shared_strains |>
    dplyr::left_join(
      meta_lookup |>
        dplyr::rename(
          subject_1 = "subject_id", pair_1 = "pair_id",
          role_1 = "role", timepoint_1 = "timepoint_days"
        ),
      by = c("sample_id_1" = "sample_id")
    ) |>
    dplyr::left_join(
      meta_lookup |>
        dplyr::rename(
          subject_2 = "subject_id", pair_2 = "pair_id",
          role_2 = "role", timepoint_2 = "timepoint_days"
        ),
      by = c("sample_id_2" = "sample_id")
    )

  # Filter to mother-child comparisons within the same pair
  mc_comparisons <- annotated |>
    dplyr::filter(
      (.data$role_1 == "mother" & .data$role_2 == "child") |
        (.data$role_1 == "child" & .data$role_2 == "mother")
    ) |>
    dplyr::filter(
      .data$pair_1 == .data$pair_2 |
        is.na(.data$pair_1) | is.na(.data$pair_2)
    )

  if (nrow(mc_comparisons) == 0L) {
    rlang::warn(
      c("No mother-child comparisons found in shared_strains.",
        "i" = "Cannot compute transmission odds."),
      class = "strainflow_no_mc_comparisons"
    )
    return(tibble::tibble(
      genome = character(),
      pair_id = character(),
      sample_id_mother = character(),
      sample_id_child = character(),
      timepoint_child = double(),
      observed_sharing = logical(),
      transmission_score = double(),
      background_rate = double(),
      prior = double(),
      transmission_odds = double(),
      log_odds = double()
    ))
  }

  # Normalize so mother is always sample_1
  mc_comparisons <- mc_comparisons |>
    dplyr::mutate(
      sample_id_mother = dplyr::if_else(
        .data$role_1 == "mother", .data$sample_id_1, .data$sample_id_2
      ),
      sample_id_child = dplyr::if_else(
        .data$role_1 == "child", .data$sample_id_1, .data$sample_id_2
      ),
      timepoint_child = dplyr::if_else(
        .data$role_1 == "child", .data$timepoint_1, .data$timepoint_2
      ),
      pair_id = dplyr::coalesce(.data$pair_1, .data$pair_2)
    )

  # -- Join with background prevalence ----------------------------------------
  mc_comparisons <- mc_comparisons |>
    dplyr::left_join(
      prevalence |> dplyr::select("genome", "background_sharing_rate"),
      by = "genome"
    ) |>
    dplyr::mutate(
      background_rate = dplyr::coalesce(.data$background_sharing_rate, 0)
    )

  # -- Compute temporal prior -------------------------------------------------
  # Prior at birth = prior_birth (0.8 by default), decays exponentially
  prior_birth <- 0.8

  mc_comparisons <- mc_comparisons |>
    dplyr::mutate(
      prior = dplyr::case_when(
        is.na(.data$timepoint_child) ~ .env$prior_weight,
        .data$timepoint_child <= 0   ~ .env$prior_birth,
        TRUE ~ .env$prior_birth * exp(-.data$timepoint_child / 365)
      ),
      # Clamp prior to avoid extremes
      prior = pmax(pmin(.data$prior, 0.99), 0.01)
    )

  # -- Compute transmission odds ----------------------------------------------
  # P(sharing | transmission) ~ transmission_score
  # P(sharing | chance) ~ background_sharing_rate
  # Bayes factor * prior odds = posterior odds

  # Use small epsilon to avoid division by zero
  epsilon <- 1e-10

  mc_comparisons <- mc_comparisons |>
    dplyr::mutate(
      # Likelihood ratio (Bayes factor)
      p_sharing_given_transmission = dplyr::coalesce(.data$transmission_score, 0),
      p_sharing_given_chance = pmax(.data$background_rate, .env$epsilon),

      # Prior odds of transmission
      prior_odds = .data$prior / (1 - .data$prior),

      # Posterior odds = Bayes factor * prior odds
      # Bayes factor = P(data|transmission) / P(data|chance)
      transmission_odds = dplyr::case_when(
        is.na(.data$transmission_score) ~ NA_real_,
        TRUE ~ (.data$p_sharing_given_transmission / .data$p_sharing_given_chance) *
          .data$prior_odds
      ),

      log_odds = dplyr::case_when(
        is.na(.data$transmission_odds) ~ NA_real_,
        .data$transmission_odds <= 0 ~ -Inf,
        TRUE ~ log10(.data$transmission_odds)
      )
    )

  # -- Format output ----------------------------------------------------------
  result <- mc_comparisons |>
    dplyr::select(
      "genome", "pair_id",
      "sample_id_mother", "sample_id_child", "timepoint_child",
      observed_sharing = "is_shared",
      "transmission_score", "background_rate",
      "prior", "transmission_odds", "log_odds"
    )

  n_high <- sum(result$log_odds > 2, na.rm = TRUE)
  n_plausible <- sum(result$log_odds > 0 & result$log_odds <= 2, na.rm = TRUE)

  cli::cli_inform(c(
    "v" = "Transmission odds computed for {nrow(result)} mother-child comparisons.",
    "i" = glue::glue(
      "High confidence (log_odds > 2): {n_high}, ",
      "Plausible (0 < log_odds <= 2): {n_plausible}"
    )
  ))

  result
}

# --------------------------------------------------------------------------
# fit_transmission_model
# --------------------------------------------------------------------------

#' Fit mixed-effects transmission model
#'
#' Fits a generalized linear mixed-effects model (GLMM) to evaluate the
#' effects of clinical covariates (delivery mode, feeding mode) on strain
#' sharing probability, controlling for species-specific background sharing
#' rates and temporal effects. Uses random intercepts for mother-child pair
#' and genome to account for non-independence.
#'
#' @param shared_strains Tibble with `is_shared` for all mother-child
#'   comparisons. Must include: `genome`, `sample_id_1`, `sample_id_2`,
#'   `is_shared`.
#' @param sample_metadata Sample metadata with: `sample_id`, `subject_id`,
#'   `pair_id`, `role`, `timepoint_days`.
#' @param pair_metadata Pair metadata with covariates: `pair_id`,
#'   `delivery_mode`, `feeding_mode`.
#' @param prevalence Background prevalence from [estimate_strain_prevalence()].
#'
#' @return A named list with components:
#'   \describe{
#'     \item{model}{The fitted `glmerMod` object from [lme4::glmer()], or
#'       `NULL` if fitting failed.}
#'     \item{summary_table}{A [tibble::tibble] of fixed effect estimates
#'       with columns: `term`, `estimate`, `std_error`, `z_value`,
#'       `p_value`.}
#'     \item{random_effects}{A [tibble::tibble] summarizing random effect
#'       variances.}
#'     \item{convergence_ok}{Logical indicating whether the model converged
#'       without warnings.}
#'     \item{n_observations}{Number of observations in the model.}
#'   }
#'
#' @export
fit_transmission_model <- function(shared_strains,
                                   sample_metadata,
                                   pair_metadata,
                                   prevalence) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(shared_strains,
                 c("genome", "sample_id_1", "sample_id_2", "is_shared"),
                 "shared_strains")
  .validate_cols(sample_metadata,
                 c("sample_id", "subject_id", "pair_id", "role", "timepoint_days"),
                 "sample_metadata")
  .validate_cols(pair_metadata, c("pair_id"), "pair_metadata")
  .validate_cols(prevalence, c("genome", "background_sharing_rate"), "prevalence")

  rlang::check_installed("lme4", reason = "to fit mixed-effects transmission model")

  # -- Build model data -------------------------------------------------------
  meta_lookup <- sample_metadata |>
    dplyr::select("sample_id", "pair_id", "role", "timepoint_days")

  # Annotate with metadata
  model_data <- shared_strains |>
    dplyr::left_join(
      meta_lookup |>
        dplyr::rename(pair_1 = "pair_id", role_1 = "role", tp_1 = "timepoint_days"),
      by = c("sample_id_1" = "sample_id")
    ) |>
    dplyr::left_join(
      meta_lookup |>
        dplyr::rename(pair_2 = "pair_id", role_2 = "role", tp_2 = "timepoint_days"),
      by = c("sample_id_2" = "sample_id")
    )

  # Filter to mother-child comparisons
  model_data <- model_data |>
    dplyr::filter(
      (.data$role_1 == "mother" & .data$role_2 == "child") |
        (.data$role_1 == "child" & .data$role_2 == "mother")
    )

  if (nrow(model_data) == 0L) {
    rlang::warn(
      "No mother-child comparisons found. Cannot fit transmission model.",
      class = "strainflow_no_mc_data"
    )
    return(list(
      model = NULL,
      summary_table = tibble::tibble(
        term = character(), estimate = double(),
        std_error = double(), z_value = double(), p_value = double()
      ),
      random_effects = tibble::tibble(
        group = character(), variance = double(), std_dev = double()
      ),
      convergence_ok = FALSE,
      n_observations = 0L
    ))
  }

  # Normalize: extract pair_id and child timepoint
  model_data <- model_data |>
    dplyr::mutate(
      pair_id = dplyr::coalesce(.data$pair_1, .data$pair_2),
      timepoint_days = dplyr::if_else(
        .data$role_1 == "child", .data$tp_1, .data$tp_2
      )
    )

  # Join with pair metadata for covariates
  # Gracefully handle missing covariates
  pair_covariates <- pair_metadata |>
    dplyr::select(
      "pair_id",
      dplyr::any_of(c("delivery_mode", "feeding_mode"))
    )

  model_data <- model_data |>
    dplyr::left_join(pair_covariates, by = "pair_id")

  # Join with prevalence
  model_data <- model_data |>
    dplyr::left_join(
      prevalence |> dplyr::select("genome", "background_sharing_rate"),
      by = "genome"
    ) |>
    dplyr::mutate(
      background_sharing_rate = dplyr::coalesce(.data$background_sharing_rate, 0),
      log_bg_rate = log(.data$background_sharing_rate + 1e-6),
      scaled_timepoint = as.numeric(scale(.data$timepoint_days)),
      is_shared_numeric = as.integer(.data$is_shared)
    )

  # Remove rows with missing outcome
  model_data <- model_data |>
    dplyr::filter(!is.na(.data$is_shared_numeric))

  # -- Check for sufficient variation -----------------------------------------
  n_shared <- sum(model_data$is_shared_numeric == 1L)
  n_not_shared <- sum(model_data$is_shared_numeric == 0L)

  if (n_shared < 5L || n_not_shared < 5L) {
    rlang::warn(
      c(
        glue::glue(
          "Insufficient variation in outcome for GLMM (shared: {n_shared}, ",
          "not shared: {n_not_shared})."
        ),
        "i" = "Need at least 5 events in each category. Returning NULL model."
      ),
      class = "strainflow_insufficient_variation"
    )
    return(list(
      model = NULL,
      summary_table = tibble::tibble(
        term = character(), estimate = double(),
        std_error = double(), z_value = double(), p_value = double()
      ),
      random_effects = tibble::tibble(
        group = character(), variance = double(), std_dev = double()
      ),
      convergence_ok = FALSE,
      n_observations = nrow(model_data)
    ))
  }

  # -- Build formula dynamically based on available covariates ----------------
  fixed_terms <- "log_bg_rate + scaled_timepoint"
  if ("delivery_mode" %in% colnames(model_data) &&
      length(unique(stats::na.omit(model_data$delivery_mode))) > 1L) {
    fixed_terms <- paste0(fixed_terms, " + delivery_mode")
  }
  if ("feeding_mode" %in% colnames(model_data) &&
      length(unique(stats::na.omit(model_data$feeding_mode))) > 1L) {
    fixed_terms <- paste0(fixed_terms, " + feeding_mode")
  }

  # Determine random effects: need >= 2 levels per grouping
  random_terms <- character()
  if (length(unique(model_data$pair_id)) >= 2L) {
    random_terms <- c(random_terms, "(1|pair_id)")
  }
  if (length(unique(model_data$genome)) >= 2L) {
    random_terms <- c(random_terms, "(1|genome)")
  }

  if (length(random_terms) == 0L) {
    # Fall back to fixed-effects only logistic regression
    formula_str <- paste0("is_shared_numeric ~ ", fixed_terms)
    rlang::inform(
      c("!" = "Insufficient grouping levels for random effects.",
        "i" = "Fitting fixed-effects logistic regression instead."),
      class = "strainflow_fixed_effects_fallback"
    )
    use_glmer <- FALSE
  } else {
    formula_str <- paste0(
      "is_shared_numeric ~ ", fixed_terms, " + ",
      paste(random_terms, collapse = " + ")
    )
    use_glmer <- TRUE
  }

  model_formula <- stats::as.formula(formula_str)

  # -- Fit model with convergence handling ------------------------------------
  convergence_ok <- TRUE
  model <- NULL

  fit_result <- tryCatch(
    {
      if (use_glmer) {
        withCallingHandlers(
          lme4::glmer(
            formula = model_formula,
            data = model_data,
            family = stats::binomial(),
            control = lme4::glmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 100000)
            )
          ),
          warning = function(w) {
            if (grepl("convergence|singular", conditionMessage(w), ignore.case = TRUE)) {
              convergence_ok <<- FALSE
              rlang::inform(
                c("!" = glue::glue("Model convergence warning: {conditionMessage(w)}")),
                class = "strainflow_convergence_warning"
              )
              invokeRestart("muffleWarning")
            }
          }
        )
      } else {
        stats::glm(
          formula = model_formula,
          data = model_data,
          family = stats::binomial()
        )
      }
    },
    error = function(e) {
      rlang::warn(
        c(
          "Transmission model fitting failed.",
          "x" = conditionMessage(e),
          "i" = "Returning NULL model. Consider simplifying covariates."
        ),
        class = "strainflow_model_fit_error"
      )
      NULL
    }
  )

  model <- fit_result

  # -- Extract results --------------------------------------------------------
  if (is.null(model)) {
    return(list(
      model = NULL,
      summary_table = tibble::tibble(
        term = character(), estimate = double(),
        std_error = double(), z_value = double(), p_value = double()
      ),
      random_effects = tibble::tibble(
        group = character(), variance = double(), std_dev = double()
      ),
      convergence_ok = FALSE,
      n_observations = nrow(model_data)
    ))
  }

  # Extract fixed effects
  summary_table <- tryCatch(
    {
      coef_summary <- summary(model)$coefficients
      tibble::tibble(
        term = rownames(coef_summary),
        estimate = coef_summary[, "Estimate"],
        std_error = coef_summary[, "Std. Error"],
        z_value = coef_summary[, if ("z value" %in% colnames(coef_summary)) "z value" else "t value"],
        p_value = coef_summary[, if ("Pr(>|z|)" %in% colnames(coef_summary)) "Pr(>|z|)" else "Pr(>|t|)"]
      )
    },
    error = function(e) {
      rlang::warn(
        c("Failed to extract model coefficients.", "x" = conditionMessage(e)),
        class = "strainflow_coef_extraction_error"
      )
      tibble::tibble(
        term = character(), estimate = double(),
        std_error = double(), z_value = double(), p_value = double()
      )
    }
  )

  # Extract random effects
  random_effects <- tryCatch(
    {
      if (use_glmer) {
        vc <- lme4::VarCorr(model)
        vc_df <- as.data.frame(vc)
        tibble::tibble(
          group = vc_df$grp,
          variance = vc_df$vcov,
          std_dev = vc_df$sdcor
        )
      } else {
        tibble::tibble(
          group = "none (fixed-effects only)",
          variance = NA_real_,
          std_dev = NA_real_
        )
      }
    },
    error = function(e) {
      tibble::tibble(
        group = character(), variance = double(), std_dev = double()
      )
    }
  )

  cli::cli_inform(c(
    "v" = "Transmission model fitted successfully ({nrow(model_data)} observations).",
    "i" = glue::glue(
      "Fixed effects: {nrow(summary_table)} terms. ",
      "Convergence: {if (convergence_ok) 'OK' else 'warnings detected'}."
    )
  ))

  list(
    model = model,
    summary_table = summary_table,
    random_effects = random_effects,
    convergence_ok = convergence_ok,
    n_observations = nrow(model_data)
  )
}

# --------------------------------------------------------------------------
# classify_transmission_event
# --------------------------------------------------------------------------

#' Classify transmission events by confidence
#'
#' Assigns a confidence classification to each potential transmission event
#' based on the Bayesian log-odds and optionally the mixed-effects model
#' results. Generates a human-readable evidence summary for each event.
#'
#' @param transmission_odds Tibble from [compute_transmission_odds()].
#' @param model_results List from [fit_transmission_model()], or `NULL`
#'   if no model was fitted.
#'
#' @return The input tibble with additional columns:
#'   \describe{
#'     \item{transmission_class}{One of `"high_confidence"`, `"plausible"`,
#'       or `"unlikely"`.}
#'     \item{evidence_summary}{Character string describing the evidence.}
#'   }
#'
#' @export
classify_transmission_event <- function(transmission_odds,
                                        model_results = NULL) {

  # -- Validate inputs --------------------------------------------------------
  .validate_cols(transmission_odds,
                 c("genome", "pair_id", "log_odds", "transmission_score",
                   "background_rate", "timepoint_child"),
                 "transmission_odds")

  # -- Classify by log-odds ---------------------------------------------------
  result <- transmission_odds |>
    dplyr::mutate(
      transmission_class = dplyr::case_when(
        is.na(.data$log_odds) ~ NA_character_,
        .data$log_odds > 2   ~ "high_confidence",
        .data$log_odds > 0   ~ "plausible",
        TRUE                 ~ "unlikely"
      )
    )

  # -- Build evidence summary -------------------------------------------------
  result <- result |>
    dplyr::mutate(
      evidence_summary = purrr::pmap_chr(
        list(
          genome = .data$genome,
          transmission_class = .data$transmission_class,
          log_odds = .data$log_odds,
          transmission_score = .data$transmission_score,
          background_rate = .data$background_rate,
          timepoint_child = .data$timepoint_child,
          observed_sharing = .data$observed_sharing
        ),
        function(genome, transmission_class, log_odds, transmission_score,
                 background_rate, timepoint_child, observed_sharing) {

          if (is.na(transmission_class)) {
            return("Insufficient data for classification.")
          }

          # Format odds ratio
          odds_str <- if (is.finite(log_odds)) {
            format(round(10^log_odds, 1), big.mark = ",")
          } else if (log_odds == Inf) {
            ">1,000,000"
          } else {
            "<0.001"
          }

          # Timepoint description
          tp_str <- if (is.na(timepoint_child)) {
            "unknown timepoint"
          } else if (timepoint_child <= 7) {
            glue::glue("birth/early neonatal (day {timepoint_child})")
          } else if (timepoint_child <= 30) {
            glue::glue("neonatal (day {timepoint_child})")
          } else if (timepoint_child <= 365) {
            glue::glue("{round(timepoint_child / 30, 1)} months")
          } else {
            glue::glue("{round(timepoint_child / 365, 1)} years")
          }

          bg_str <- format(round(background_rate * 100, 2), nsmall = 2)
          score_str <- format(round(transmission_score, 6), nsmall = 6)

          sharing_str <- if (isTRUE(observed_sharing)) {
            "Strain sharing detected."
          } else {
            "No strain sharing detected."
          }

          glue::glue(
            "{sharing_str} ",
            "{genome} at {tp_str}: ",
            "transmission score={score_str}, ",
            "background prevalence={bg_str}%, ",
            "odds ratio={odds_str}x favoring transmission. ",
            "Classification: {transmission_class}."
          )
        }
      )
    )

  # -- Summarize if model results available -----------------------------------
  if (!is.null(model_results) && !is.null(model_results$model)) {
    model_note <- tryCatch(
      {
        st <- model_results$summary_table
        sig_terms <- st |>
          dplyr::filter(.data$p_value < 0.05, .data$term != "(Intercept)")
        if (nrow(sig_terms) > 0L) {
          glue::glue(
            " GLMM significant effects: {paste(sig_terms$term, collapse = ', ')}."
          )
        } else {
          " GLMM: no significant fixed effects (p < 0.05)."
        }
      },
      error = function(e) ""
    )

    if (nzchar(model_note)) {
      result <- result |>
        dplyr::mutate(
          evidence_summary = paste0(.data$evidence_summary, .env$model_note)
        )
    }
  }

  n_high <- sum(result$transmission_class == "high_confidence", na.rm = TRUE)
  n_plausible <- sum(result$transmission_class == "plausible", na.rm = TRUE)
  n_unlikely <- sum(result$transmission_class == "unlikely", na.rm = TRUE)

  cli::cli_inform(c(
    "v" = "Transmission events classified: {nrow(result)} events.",
    "i" = glue::glue(
      "high_confidence: {n_high}, plausible: {n_plausible}, unlikely: {n_unlikely}"
    )
  ))

  result
}
