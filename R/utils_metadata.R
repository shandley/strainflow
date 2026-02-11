# ============================================================================
# utils_metadata.R - Sample and pair metadata loading, validation, and
#                    comparison manifest generation
# ============================================================================

#' Required columns in sample metadata
#' @noRd
.sample_required_cols <- c(
  "sample_id", "subject_id", "role", "timepoint_days",
  "pair_id", "fastq_r1", "fastq_r2"
)

#' Optional columns in sample metadata
#' @noRd
.sample_optional_cols <- c("delivery_mode", "feeding_mode", "antibiotics")

#' Allowed values for the role column
#' @noRd
.valid_roles <- c("mother", "child")

#' Required columns in pair metadata
#' @noRd
.pair_required_cols <- c("pair_id", "mother_subject_id", "child_subject_id")

# --------------------------------------------------------------------------
# load_sample_metadata
# --------------------------------------------------------------------------

#' Load and validate sample metadata
#'
#' Reads a CSV sample sheet and validates column names, data types, uniqueness,
#' allowed factor levels, and FASTQ file existence.
#'
#' @param sample_sheet_path Path to \code{sample_metadata.csv}.
#' @return A validated \code{\link[tibble]{tibble}} sorted by \code{pair_id},
#'   \code{role}, and \code{timepoint_days}.
#' @export
load_sample_metadata <- function(sample_sheet_path) {
  if (!fs::file_exists(sample_sheet_path)) {
    rlang::abort(
      glue::glue("Sample metadata file not found: {sample_sheet_path}"),
      class = "strainflow_metadata_missing"
    )
  }

  sm <- readr::read_csv(
    sample_sheet_path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )

  # --- required columns ---
  present <- colnames(sm)
  missing <- setdiff(.sample_required_cols, present)
  if (length(missing) > 0L) {
    rlang::abort(
      c(
        glue::glue(
          "Sample metadata is missing required columns: ",
          "{paste(missing, collapse = ', ')}"
        ),
        "i" = glue::glue(
          "Required: {paste(.sample_required_cols, collapse = ', ')}"
        )
      ),
      class = "strainflow_metadata_cols"
    )
  }

  # --- coerce types ---
  sm <- sm |>
    dplyr::mutate(
      timepoint_days = suppressWarnings(as.numeric(.data$timepoint_days))
    )

  # --- validate unique sample_id ---
  dup_ids <- sm$sample_id[duplicated(sm$sample_id)]
  if (length(dup_ids) > 0L) {
    rlang::abort(
      glue::glue(
        "Duplicate sample_id values found: {paste(unique(dup_ids), collapse = ', ')}"
      ),
      class = "strainflow_metadata_dup_ids"
    )
  }

  # --- validate role ---
  invalid_roles <- setdiff(unique(sm$role), .valid_roles)
  if (length(invalid_roles) > 0L) {
    rlang::abort(
      glue::glue(
        "Invalid role values: {paste(invalid_roles, collapse = ', ')}. ",
        "Allowed values: {paste(.valid_roles, collapse = ', ')}"
      ),
      class = "strainflow_metadata_roles"
    )
  }

  # --- validate timepoint_days is numeric ---
  if (any(is.na(sm$timepoint_days))) {
    bad_rows <- which(is.na(sm$timepoint_days))
    rlang::abort(
      glue::glue(
        "Non-numeric or missing timepoint_days at rows: ",
        "{paste(bad_rows, collapse = ', ')}"
      ),
      class = "strainflow_metadata_timepoints"
    )
  }

  # --- validate FASTQ paths exist ---
  fastq_paths <- c(sm$fastq_r1, sm$fastq_r2)
  missing_fq <- fastq_paths[!fs::file_exists(fastq_paths)]
  if (length(missing_fq) > 0L) {
    # Show at most 10 paths to avoid flooding the console
    shown <- utils::head(missing_fq, 10L)
    suffix <- if (length(missing_fq) > 10L) {
      glue::glue(" ... and {length(missing_fq) - 10L} more")
    } else {
      ""
    }
    rlang::abort(
      c(
        glue::glue(
          "{length(missing_fq)} FASTQ file(s) not found."
        ),
        "x" = paste0(shown, collapse = "\n"),
        "i" = suffix
      ),
      class = "strainflow_metadata_fastq_missing"
    )
  }

  # --- sort and return ---
  sm <- sm |>
    dplyr::arrange(.data$pair_id, .data$role, .data$timepoint_days)

  cli::cli_inform(c(
    "v" = "Sample metadata loaded: {nrow(sm)} samples across {length(unique(sm$pair_id))} pairs."
  ))

  sm
}

# --------------------------------------------------------------------------
# load_pair_metadata
# --------------------------------------------------------------------------

#' Load pair metadata
#'
#' Reads and validates a pair-level metadata CSV that links mother-child pairs.
#'
#' @param pair_sheet_path Path to \code{pair_metadata.csv}.
#' @return A validated \code{\link[tibble]{tibble}} with at least
#'   \code{pair_id}, \code{mother_subject_id}, and \code{child_subject_id}.
#' @export
load_pair_metadata <- function(pair_sheet_path) {
  if (!fs::file_exists(pair_sheet_path)) {
    rlang::abort(
      glue::glue("Pair metadata file not found: {pair_sheet_path}"),
      class = "strainflow_metadata_missing"
    )
  }

  pm <- readr::read_csv(
    pair_sheet_path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )

  missing <- setdiff(.pair_required_cols, colnames(pm))
  if (length(missing) > 0L) {
    rlang::abort(
      c(
        glue::glue(
          "Pair metadata is missing required columns: ",
          "{paste(missing, collapse = ', ')}"
        ),
        "i" = glue::glue(
          "Required: {paste(.pair_required_cols, collapse = ', ')}"
        )
      ),
      class = "strainflow_metadata_cols"
    )
  }

  # --- validate unique pair_id ---
  dup_ids <- pm$pair_id[duplicated(pm$pair_id)]
  if (length(dup_ids) > 0L) {
    rlang::abort(
      glue::glue(
        "Duplicate pair_id values in pair metadata: ",
        "{paste(unique(dup_ids), collapse = ', ')}"
      ),
      class = "strainflow_metadata_dup_ids"
    )
  }

  cli::cli_inform(c(
    "v" = "Pair metadata loaded: {nrow(pm)} mother-child pairs."
  ))

  pm
}

# --------------------------------------------------------------------------
# generate_comparison_manifest
# --------------------------------------------------------------------------

#' Generate comparison manifest for inStrain compare
#'
#' Builds a table of all pairwise sample comparisons needed for inStrain
#' compare, including:
#' \itemize{
#'   \item \strong{mother_child}: every mother sample vs. every child sample
#'     within the same pair (all timepoint combinations).
#'   \item \strong{within_mother}: longitudinal comparisons of all mother
#'     samples within a pair.
#'   \item \strong{within_child}: longitudinal comparisons of all child
#'     samples within a pair.
#' }
#'
#' @param sample_metadata A tibble from \code{\link{load_sample_metadata}}.
#' @param pair_metadata A tibble from \code{\link{load_pair_metadata}}.
#' @return A \code{\link[tibble]{tibble}} with columns \code{comparison_id},
#'   \code{sample_id_1}, \code{sample_id_2}, \code{pair_id}, and
#'   \code{comparison_type}.
#' @export
generate_comparison_manifest <- function(sample_metadata, pair_metadata) {
  all_comparisons <- purrr::map_dfr(
    unique(pair_metadata$pair_id),
    function(pid) {
      pair_samples <- sample_metadata |>
        dplyr::filter(.data$pair_id == .env$pid)

      mothers <- pair_samples |>
        dplyr::filter(.data$role == "mother") |>
        dplyr::pull(.data$sample_id)

      children <- pair_samples |>
        dplyr::filter(.data$role == "child") |>
        dplyr::pull(.data$sample_id)

      # -- mother vs child (all combinations) --
      mc_comp <- if (length(mothers) > 0L && length(children) > 0L) {
        tidyr::expand_grid(
          sample_id_1 = mothers,
          sample_id_2 = children
        ) |>
          dplyr::mutate(
            pair_id = pid,
            comparison_type = "mother_child"
          )
      } else {
        tibble::tibble(
          sample_id_1 = character(),
          sample_id_2 = character(),
          pair_id = character(),
          comparison_type = character()
        )
      }

      # -- within-mother longitudinal --
      wm_comp <- .generate_within_comparisons(mothers, pid, "within_mother")

      # -- within-child longitudinal --
      wc_comp <- .generate_within_comparisons(children, pid, "within_child")

      dplyr::bind_rows(mc_comp, wm_comp, wc_comp)
    }
  )

  # Build readable comparison_id
  all_comparisons <- all_comparisons |>
    dplyr::mutate(
      comparison_id = glue::glue(
        "{pair_id}__{comparison_type}__{sample_id_1}_vs_{sample_id_2}"
      )
    ) |>
    dplyr::select(
      "comparison_id", "sample_id_1", "sample_id_2",
      "pair_id", "comparison_type"
    )

  cli::cli_inform(c(
    "v" = "Comparison manifest generated: {nrow(all_comparisons)} comparisons.",
    "i" = paste0(
      "mother_child: ",
      sum(all_comparisons$comparison_type == "mother_child"),
      ", within_mother: ",
      sum(all_comparisons$comparison_type == "within_mother"),
      ", within_child: ",
      sum(all_comparisons$comparison_type == "within_child")
    )
  ))

  all_comparisons
}

#' Generate within-group pairwise comparisons (longitudinal)
#'
#' @param sample_ids Character vector of sample IDs within a group.
#' @param pid Pair ID string.
#' @param comp_type Comparison type label.
#' @return Tibble of pairwise combinations.
#' @noRd
.generate_within_comparisons <- function(sample_ids, pid, comp_type) {
  if (length(sample_ids) < 2L) {
    return(tibble::tibble(
      sample_id_1 = character(),
      sample_id_2 = character(),
      pair_id = character(),
      comparison_type = character()
    ))
  }

  pairs <- utils::combn(sort(sample_ids), 2L, simplify = FALSE)

  tibble::tibble(
    sample_id_1     = purrr::map_chr(pairs, 1L),
    sample_id_2     = purrr::map_chr(pairs, 2L),
    pair_id         = pid,
    comparison_type = comp_type
  )
}

# --------------------------------------------------------------------------
# generate_negative_controls
# --------------------------------------------------------------------------

#' Generate negative control comparisons from unrelated pairs
#'
#' For each true mother-child pair, randomly selects \code{n_controls} unrelated
#' mother-child combinations to serve as negative controls. This is critical for
#' calibrating transmission thresholds.
#'
#' @param sample_metadata A tibble from \code{\link{load_sample_metadata}}.
#' @param pair_metadata A tibble from \code{\link{load_pair_metadata}}.
#' @param n_controls Number of random unrelated pairs to generate per true pair.
#'   Defaults to 5.
#' @return A \code{\link[tibble]{tibble}} with the same structure as
#'   \code{\link{generate_comparison_manifest}} but with
#'   \code{comparison_type = "negative_control"}.
#' @export
generate_negative_controls <- function(sample_metadata,
                                       pair_metadata,
                                       n_controls = 5L) {
  n_controls <- as.integer(n_controls)
  if (n_controls < 1L) {
    rlang::abort(
      "n_controls must be a positive integer.",
      class = "strainflow_invalid_arg"
    )
  }

  pair_ids <- unique(pair_metadata$pair_id)

  if (length(pair_ids) < 2L) {
    rlang::abort(
      "Need at least 2 pairs to generate negative controls.",
      class = "strainflow_insufficient_pairs"
    )
  }

  # Build a lookup: pair_id -> mother sample_ids, child sample_ids
  pair_lookup <- purrr::map(
    stats::setNames(pair_ids, pair_ids),
    function(pid) {
      pair_samples <- sample_metadata |>
        dplyr::filter(.data$pair_id == .env$pid)
      list(
        mothers = pair_samples |>
          dplyr::filter(.data$role == "mother") |>
          dplyr::pull(.data$sample_id),
        children = pair_samples |>
          dplyr::filter(.data$role == "child") |>
          dplyr::pull(.data$sample_id)
      )
    }
  )

  neg_controls <- purrr::map_dfr(pair_ids, function(focal_pid) {
    # Get other pair_ids
    other_pids <- setdiff(pair_ids, focal_pid)

    # Sample up to n_controls unrelated partners
    n_avail <- length(other_pids)
    n_pick <- min(n_controls, n_avail)
    selected_pids <- sample(other_pids, size = n_pick)

    focal_mothers <- pair_lookup[[focal_pid]]$mothers
    # For each selected unrelated pair, cross the focal mother with unrelated child
    purrr::map_dfr(selected_pids, function(other_pid) {
      other_children <- pair_lookup[[other_pid]]$children

      if (length(focal_mothers) == 0L || length(other_children) == 0L) {
        return(tibble::tibble(
          sample_id_1 = character(),
          sample_id_2 = character(),
          pair_id = character(),
          comparison_type = character(),
          true_pair_id = character(),
          unrelated_pair_id = character()
        ))
      }

      # Pick one representative mother and child to keep the control set manageable
      mother_rep <- sample(focal_mothers, 1L)
      child_rep  <- sample(other_children, 1L)

      tibble::tibble(
        sample_id_1     = mother_rep,
        sample_id_2     = child_rep,
        pair_id         = glue::glue("{focal_pid}_vs_{other_pid}"),
        comparison_type = "negative_control",
        true_pair_id    = focal_pid,
        unrelated_pair_id = other_pid
      )
    })
  })

  # Build comparison_id
  neg_controls <- neg_controls |>
    dplyr::mutate(
      comparison_id = glue::glue(
        "neg__{true_pair_id}_vs_{unrelated_pair_id}__{sample_id_1}_vs_{sample_id_2}"
      )
    ) |>
    dplyr::select(
      "comparison_id", "sample_id_1", "sample_id_2",
      "pair_id", "comparison_type",
      "true_pair_id", "unrelated_pair_id"
    )

  cli::cli_inform(c(
    "v" = "Negative controls generated: {nrow(neg_controls)} unrelated comparisons.",
    "i" = "{n_controls} controls requested per pair ({length(pair_ids)} pairs)."
  ))

  neg_controls
}

# --------------------------------------------------------------------------
# build_pair_profile_groups
# --------------------------------------------------------------------------

#' Build profile groupings for inStrain compare by pair
#'
#' Creates a table mapping each mother-child pair to the list of inStrain
#' profile directories for all samples belonging to that pair. This is used
#' to assemble the input for \code{inStrain compare}.
#'
#' @param instrain_results A tibble or data frame with at least columns
#'   \code{sample_id} and \code{profile_dir} (path to the inStrain profile
#'   output directory for each sample).
#' @param sample_metadata A tibble from \code{\link{load_sample_metadata}}.
#' @param pair_metadata A tibble from \code{\link{load_pair_metadata}}.
#' @return A \code{\link[tibble]{tibble}} with columns \code{pair_id} and
#'   \code{profile_dirs} (a list-column of character vectors of profile
#'   directory paths).
#' @export
build_pair_profile_groups <- function(instrain_results,
                                      sample_metadata,
                                      pair_metadata) {
  # Validate instrain_results columns
  required_cols <- c("sample_id", "profile_dir")
  missing_cols <- setdiff(required_cols, colnames(instrain_results))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "instrain_results is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_results_cols"
    )
  }

  # Validate profile directories exist
  missing_dirs <- instrain_results$profile_dir[
    !fs::dir_exists(instrain_results$profile_dir)
  ]
  if (length(missing_dirs) > 0L) {
    shown <- utils::head(missing_dirs, 5L)
    rlang::warn(
      c(
        glue::glue(
          "{length(missing_dirs)} inStrain profile director{ifelse(length(missing_dirs) == 1, 'y does', 'ies do')} not exist."
        ),
        "!" = paste0(shown, collapse = "\n")
      ),
      class = "strainflow_profiles_missing"
    )
  }

  # Join sample metadata to get pair_id, then group
  profile_with_pair <- instrain_results |>
    dplyr::select("sample_id", "profile_dir") |>
    dplyr::inner_join(
      sample_metadata |> dplyr::select("sample_id", "pair_id"),
      by = "sample_id"
    )

  # Filter to only pairs present in pair_metadata
  valid_pairs <- unique(pair_metadata$pair_id)
  profile_with_pair <- profile_with_pair |>
    dplyr::filter(.data$pair_id %in% .env$valid_pairs)

  # Group into list-column
  grouped <- profile_with_pair |>
    dplyr::group_by(.data$pair_id) |>
    dplyr::summarise(
      profile_dirs = list(.data$profile_dir),
      n_profiles = dplyr::n(),
      .groups = "drop"
    )

  cli::cli_inform(c(
    "v" = "Profile groups built for {nrow(grouped)} pairs ",
    "i" = "Total profiles: {sum(grouped$n_profiles)}."
  ))

  grouped
}
