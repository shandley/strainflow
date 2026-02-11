# ============================================================================
# step11_visualization.R - Publication-quality visualization functions
# ============================================================================
#
# All 12 plotting functions plus a package theme helper. Each function takes
# data tibbles as input and returns a ggplot object (or saves to file and
# returns the path). Uses ggplot2 as the primary plotting library with a
# consistent theme across all figures.
# ============================================================================

# --------------------------------------------------------------------------
# Consistent color palettes used across the package
# --------------------------------------------------------------------------

#' Strainflow color palettes (internal)
#' @noRd
.sf_colors <- list(
  related       = "#D73027",
  unrelated     = "#969696",
  mother        = "#4575B4",
  child         = "#FC8D59",
  host          = "#D73027",

  nonhost       = "#4575B4",
  reads_before  = "#BDBDBD",
  reads_after   = "#4575B4",
  significant   = "#D73027",
  nonsignificant = "#969696",
  hypermutator       = "#E66101",
  hyper_recombinator = "#5E3C99",
  balanced           = "#969696",
  shared_strain = "#D73027",
  not_shared    = "#F0F0F0"
)


# --------------------------------------------------------------------------
# theme_strainflow
# --------------------------------------------------------------------------

#' StrainFlow default ggplot theme
#'
#' Provides a clean, publication-ready theme used consistently across all
#' strainflow visualization functions.
#'
#' @return A \code{\link[ggplot2]{theme}} object.
#' @export
theme_strainflow <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40"),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 9)
    )
}


# --------------------------------------------------------------------------
# Internal helper: conditionally save or return
# --------------------------------------------------------------------------

#' Save a ggplot to file or return the object
#'
#' @param p A ggplot object.
#' @param output_path Optional file path. If provided, saves and returns the
#'   path; otherwise returns the ggplot.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @return ggplot object or character scalar file path.
#' @noRd
.save_or_return <- function(p, output_path = NULL, width = 10, height = 7) {
  if (!is.null(output_path)) {
    fs::dir_create(fs::path_dir(output_path), recurse = TRUE)
    ggplot2::ggsave(
      filename = output_path,
      plot = p,
      width = width,
      height = height,
      dpi = 300
    )
    return(as.character(output_path))
  }
  p
}


# ==========================================================================
# QC Plots
# ==========================================================================

# --------------------------------------------------------------------------
# 1. plot_qc_summary
# --------------------------------------------------------------------------

#' Plot QC summary: reads before and after filtering
#'
#' Produces a grouped bar chart showing the number of reads before (grey) and
#' after (blue) quality filtering for each sample. Optionally facets by
#' mother/child role when sample metadata is provided.
#'
#' @param qc_stats A \code{\link[tibble]{tibble}} with columns \code{sample_id},
#'   \code{reads_before}, \code{reads_after}, and \code{pct_passing}.
#' @param sample_metadata Optional sample metadata tibble with \code{sample_id},
#'   \code{role}, \code{pair_id}, and \code{timepoint_days} for ordering and
#'   faceting.
#' @param output_path Optional file path to save the plot. If \code{NULL}, the
#'   ggplot object is returned.
#' @return A ggplot object or the file path (character scalar) if saved.
#' @export
plot_qc_summary <- function(qc_stats,
                            sample_metadata = NULL,
                            output_path = NULL) {

  # --- validate required columns ---
  required <- c("sample_id", "reads_before", "reads_after")
  missing_cols <- setdiff(required, colnames(qc_stats))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "qc_stats is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  # --- optionally merge metadata for ordering and faceting ---
  plot_data <- qc_stats
  has_meta <- !is.null(sample_metadata) &&
    all(c("sample_id", "role", "pair_id", "timepoint_days") %in%
        colnames(sample_metadata))

  if (has_meta) {
    plot_data <- plot_data |>
      dplyr::left_join(
        sample_metadata |>
          dplyr::select("sample_id", "role", "pair_id", "timepoint_days"),
        by = "sample_id"
      ) |>
      dplyr::arrange(.data$pair_id, .data$timepoint_days)
  }

  # --- pivot to long format for grouped bars ---
  plot_long <- plot_data |>
    tidyr::pivot_longer(
      cols = c("reads_before", "reads_after"),
      names_to = "stage",
      values_to = "read_count"
    ) |>
    dplyr::mutate(
      stage = factor(
        .data$stage,
        levels = c("reads_before", "reads_after"),
        labels = c("Before QC", "After QC")
      )
    )

  # Preserve sample ordering in the factor
  sample_order <- unique(plot_data$sample_id)
  plot_long <- plot_long |>
    dplyr::mutate(
      sample_id = factor(.data$sample_id, levels = sample_order)
    )

  # --- compute median after-QC reads for reference line ---
  median_after <- stats::median(qc_stats$reads_after, na.rm = TRUE)

  # --- build plot ---
  p <- ggplot2::ggplot(
    plot_long,
    ggplot2::aes(x = .data$sample_id, y = .data$read_count, fill = .data$stage)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = 0.7) +
    ggplot2::geom_hline(
      yintercept = median_after,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.5
    ) +
    ggplot2::scale_fill_manual(
      values = c("Before QC" = .sf_colors$reads_before,
                 "After QC"  = .sf_colors$reads_after),
      name = NULL
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      title = "QC Summary: Read Counts Before and After Filtering",
      subtitle = glue::glue(
        "Dashed line = median post-QC reads ({scales::comma(median_after)})"
      ),
      x = "Sample",
      y = "Number of Reads"
    ) +
    theme_strainflow() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
    )

  # --- facet by role if metadata is available ---
  if (has_meta) {
    p <- p +
      ggplot2::facet_wrap(
        ~ .data$role,
        scales = "free_x",
        ncol = 2
      )
  }

  .save_or_return(p, output_path, width = 12, height = 6)
}


# --------------------------------------------------------------------------
# 2. plot_host_removal
# --------------------------------------------------------------------------

#' Plot host removal rates
#'
#' Produces a stacked bar chart showing host reads (red) and non-host reads
#' (blue) per sample, sorted by percentage of host reads in descending order.
#'
#' @param host_stats A \code{\link[tibble]{tibble}} with columns
#'   \code{sample_id}, \code{pct_host}, \code{nonhost_reads}, and
#'   \code{host_reads}.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_host_removal <- function(host_stats, output_path = NULL) {

  required <- c("sample_id", "pct_host", "nonhost_reads", "host_reads")
  missing_cols <- setdiff(required, colnames(host_stats))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "host_stats is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  # --- sort samples by pct_host descending ---
  plot_data <- host_stats |>
    dplyr::arrange(dplyr::desc(.data$pct_host)) |>
    dplyr::mutate(
      sample_id = factor(.data$sample_id, levels = .data$sample_id)
    )

  # --- pivot to long format ---
  plot_long <- plot_data |>
    tidyr::pivot_longer(
      cols = c("host_reads", "nonhost_reads"),
      names_to = "read_type",
      values_to = "read_count"
    ) |>
    dplyr::mutate(
      read_type = factor(
        .data$read_type,
        levels = c("host_reads", "nonhost_reads"),
        labels = c("Host", "Non-host")
      )
    )

  p <- ggplot2::ggplot(
    plot_long,
    ggplot2::aes(x = .data$sample_id, y = .data$read_count,
                 fill = .data$read_type)
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(
      values = c("Host"     = .sf_colors$host,
                 "Non-host" = .sf_colors$nonhost),
      name = NULL
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      title = "Host Removal: Read Composition per Sample",
      subtitle = "Samples sorted by percentage of host reads (descending)",
      x = "Sample",
      y = "Number of Reads"
    ) +
    theme_strainflow() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
    )

  .save_or_return(p, output_path, width = 12, height = 6)
}


# --------------------------------------------------------------------------
# 3. plot_taxonomic_composition
# --------------------------------------------------------------------------

#' Plot taxonomic composition
#'
#' Produces a stacked bar chart of relative abundances for each sample. The top
#' N species (by mean abundance across all samples) are shown individually;
#' remaining taxa are grouped as "Other".
#'
#' @param abundance_long A long-format \code{\link[tibble]{tibble}} with
#'   columns \code{sample_id}, \code{clade_name}, and
#'   \code{relative_abundance}.
#' @param sample_metadata Optional sample metadata tibble with
#'   \code{sample_id}, \code{pair_id}, \code{role}, and
#'   \code{timepoint_days} for faceting and ordering.
#' @param top_n Integer. Number of top species to display individually.
#'   Defaults to 20.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_taxonomic_composition <- function(abundance_long,
                                       sample_metadata = NULL,
                                       top_n = 20L,
                                       output_path = NULL) {

  required <- c("sample_id", "clade_name", "relative_abundance")
  missing_cols <- setdiff(required, colnames(abundance_long))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "abundance_long is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  top_n <- as.integer(top_n)

  # --- identify top N species by mean relative abundance ---
  mean_abund <- abundance_long |>
    dplyr::group_by(.data$clade_name) |>
    dplyr::summarise(
      mean_ra = mean(.data$relative_abundance, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$mean_ra))

  top_species <- utils::head(mean_abund$clade_name, top_n)

  # --- group remaining species as "Other" ---
  plot_data <- abundance_long |>
    dplyr::mutate(
      taxon = ifelse(
        .data$clade_name %in% top_species,
        .data$clade_name,
        "Other"
      )
    ) |>
    dplyr::group_by(.data$sample_id, .data$taxon) |>
    dplyr::summarise(
      relative_abundance = sum(.data$relative_abundance, na.rm = TRUE),
      .groups = "drop"
    )

  # --- order the taxon factor: top species by mean abundance, then "Other" ---
  taxon_levels <- c(top_species, "Other")
  plot_data <- plot_data |>
    dplyr::mutate(
      taxon = factor(.data$taxon, levels = rev(taxon_levels))
    )

  # --- merge metadata for ordering and faceting ---
  has_meta <- !is.null(sample_metadata) &&
    all(c("sample_id", "pair_id", "role", "timepoint_days") %in%
        colnames(sample_metadata))

  if (has_meta) {
    meta_order <- sample_metadata |>
      dplyr::select("sample_id", "pair_id", "role", "timepoint_days") |>
      dplyr::arrange(.data$pair_id, .data$role, .data$timepoint_days)

    plot_data <- plot_data |>
      dplyr::left_join(meta_order, by = "sample_id") |>
      dplyr::mutate(
        sample_id = factor(.data$sample_id, levels = unique(meta_order$sample_id))
      )
  }

  # --- choose color palette ---
  n_colors <- length(taxon_levels)
  if (n_colors <= 12) {
    palette <- scales::hue_pal()(n_colors)
  } else {
    # Use viridis-based palette (already a dependency) with "Other" in grey
    palette <- c(
      viridis::viridis(n_colors - 1, option = "D"),
      "#CCCCCC"
    )
  }
  names(palette) <- taxon_levels

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$sample_id, y = .data$relative_abundance,
                 fill = .data$taxon)
  ) +
    ggplot2::geom_col(width = 0.9, position = "stack") +
    ggplot2::scale_fill_manual(
      values = palette,
      name = "Species",
      guide = ggplot2::guide_legend(ncol = 2, reverse = TRUE)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(scale = 1),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = "Taxonomic Composition",
      subtitle = glue::glue("Top {top_n} species by mean relative abundance"),
      x = "Sample",
      y = "Relative Abundance (%)"
    ) +
    theme_strainflow() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      legend.text = ggplot2::element_text(size = 7),
      legend.key.size = ggplot2::unit(0.4, "cm")
    )

  if (has_meta) {
    p <- p +
      ggplot2::facet_wrap(
        ~ .data$pair_id,
        scales = "free_x",
        nrow = 1
      )
  }

  .save_or_return(p, output_path, width = 16, height = 8)
}


# ==========================================================================
# Dual-Signal Plots
# ==========================================================================

# --------------------------------------------------------------------------
# 4. plot_popani_vs_apss
# --------------------------------------------------------------------------

#' Plot popANI vs APSS scatter
#'
#' Scatter plot of population-level average nucleotide identity (popANI) versus
#' average pairwise synteny score (APSS). Related pairs are shown in red and
#' unrelated pairs in grey. Dashed lines mark the classification thresholds;
#' the upper-right quadrant represents shared strains.
#'
#' @param transmission_scores A \code{\link[tibble]{tibble}} with columns
#'   \code{popANI}, \code{apss}, and \code{comparison_type} (with values like
#'   "related" and "unrelated").
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_popani_vs_apss <- function(transmission_scores, output_path = NULL) {

  required <- c("popANI", "apss", "comparison_type")
  missing_cols <- setdiff(required, colnames(transmission_scores))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "transmission_scores is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  popani_threshold <- 0.99999
  apss_threshold   <- 0.95

  p <- ggplot2::ggplot(
    transmission_scores,
    ggplot2::aes(x = .data$popANI, y = .data$apss, color = .data$comparison_type)
  ) +
    # --- draw unrelated first so related points are on top ---
    ggplot2::geom_point(
      data = transmission_scores |>
        dplyr::filter(.data$comparison_type != "related"),
      alpha = 0.4, size = 1.5
    ) +
    ggplot2::geom_point(
      data = transmission_scores |>
        dplyr::filter(.data$comparison_type == "related"),
      alpha = 0.7, size = 2
    ) +
    ggplot2::geom_vline(
      xintercept = popani_threshold,
      linetype = "dashed", color = "grey30", linewidth = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = apss_threshold,
      linetype = "dashed", color = "grey30", linewidth = 0.5
    ) +
    ggplot2::scale_color_manual(
      values = c("related"   = .sf_colors$related,
                 "unrelated" = .sf_colors$unrelated),
      name = "Comparison Type",
      labels = c("related" = "Related", "unrelated" = "Unrelated")
    ) +
    ggplot2::annotate(
      "text",
      x = popani_threshold + (1 - popani_threshold) * 0.5,
      y = apss_threshold + (1 - apss_threshold) * 0.5,
      label = "Shared\nStrains",
      fontface = "italic",
      color = "grey40",
      size = 3.5
    ) +
    ggplot2::labs(
      title = "Dual-Signal Strain Comparison: popANI vs APSS",
      subtitle = glue::glue(
        "popANI threshold = {popani_threshold}; APSS threshold = {apss_threshold}"
      ),
      x = "popANI (Population Average Nucleotide Identity)",
      y = "APSS (Average Pairwise Synteny Score)"
    ) +
    theme_strainflow()

  .save_or_return(p, output_path, width = 9, height = 7)
}


# --------------------------------------------------------------------------
# 5. plot_evolutionary_modes
# --------------------------------------------------------------------------

#' Plot species evolutionary mode classification
#'
#' Scatter plot of SNV variance vs synteny variance, colored by evolutionary
#' mode classification. Top species are labeled with \code{ggrepel}.
#'
#' @param evolutionary_modes A \code{\link[tibble]{tibble}} from
#'   \code{classify_evolutionary_mode}, expected to contain \code{genome} (or
#'   \code{species}), \code{snv_variance}, \code{synteny_variance}, and
#'   \code{evolutionary_mode}.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_evolutionary_modes <- function(evolutionary_modes, output_path = NULL) {

  # Allow either "genome" or "species" as the identifier column
  id_col <- if ("genome" %in% colnames(evolutionary_modes)) {
    "genome"
  } else if ("species" %in% colnames(evolutionary_modes)) {
    "species"
  } else {
    rlang::abort(
      "evolutionary_modes must contain a 'genome' or 'species' column.",
      class = "strainflow_column_missing"
    )
  }

  required <- c("snv_variance", "synteny_variance", "evolutionary_mode")
  missing_cols <- setdiff(required, colnames(evolutionary_modes))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "evolutionary_modes is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  # --- compute medians for crosshair lines ---
  med_snv <- stats::median(evolutionary_modes$snv_variance, na.rm = TRUE)
  med_syn <- stats::median(evolutionary_modes$synteny_variance, na.rm = TRUE)

  # --- label top species (extreme SNV or synteny variance) ---
  top_n_label <- min(10L, nrow(evolutionary_modes))
  top_snv <- evolutionary_modes |>
    dplyr::arrange(dplyr::desc(.data$snv_variance)) |>
    utils::head(ceiling(top_n_label / 2))
  top_syn <- evolutionary_modes |>
    dplyr::arrange(dplyr::desc(.data$synteny_variance)) |>
    utils::head(ceiling(top_n_label / 2))
  label_ids <- unique(c(top_snv[[id_col]], top_syn[[id_col]]))

  evolutionary_modes <- evolutionary_modes |>
    dplyr::mutate(
      label = ifelse(.data[[id_col]] %in% label_ids, .data[[id_col]], NA_character_)
    )

  p <- ggplot2::ggplot(
    evolutionary_modes,
    ggplot2::aes(
      x = .data$snv_variance,
      y = .data$synteny_variance,
      color = .data$evolutionary_mode
    )
  ) +
    ggplot2::geom_vline(
      xintercept = med_snv, linetype = "dashed", color = "grey60",
      linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = med_syn, linetype = "dashed", color = "grey60",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(alpha = 0.7, size = 2.5) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$label),
      size = 3,
      max.overlaps = 15,
      segment.color = "grey50",
      segment.size = 0.3,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "hypermutator"       = .sf_colors$hypermutator,
        "hyper_recombinator" = .sf_colors$hyper_recombinator,
        "balanced"           = .sf_colors$balanced
      ),
      name = "Evolutionary Mode",
      labels = c(
        "hypermutator"       = "Hypermutator",
        "hyper_recombinator" = "Hyper-recombinator",
        "balanced"           = "Balanced"
      )
    ) +
    ggplot2::labs(
      title = "Species Evolutionary Mode Classification",
      subtitle = "Dashed crosshairs mark median SNV and synteny variance",
      x = "SNV Variance",
      y = "Synteny Variance"
    ) +
    theme_strainflow()

  .save_or_return(p, output_path, width = 10, height = 8)
}


# --------------------------------------------------------------------------
# 6. plot_strain_sharing_heatmap
# --------------------------------------------------------------------------

#' Plot strain sharing heatmap
#'
#' Tile-based heatmap (using \code{geom_tile}) showing strain sharing across
#' pairs. Rows represent genomes (species) and columns represent mother-child
#' pairs. Tiles are colored by transmission score (viridis scale) or binary
#' sharing status.
#'
#' @param shared_strains A \code{\link[tibble]{tibble}} with columns
#'   \code{genome}, \code{pair_id}, and either \code{transmission_score}
#'   (numeric) or \code{is_shared} (logical).
#' @param pair_metadata Optional pair metadata tibble for column annotation.
#'   If it contains \code{delivery_mode}, a top annotation strip is added.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_strain_sharing_heatmap <- function(shared_strains,
                                        pair_metadata = NULL,
                                        output_path = NULL) {

  required <- c("genome", "pair_id")
  missing_cols <- setdiff(required, colnames(shared_strains))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "shared_strains is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  has_score <- "transmission_score" %in% colnames(shared_strains)
  has_binary <- "is_shared" %in% colnames(shared_strains)

  if (!has_score && !has_binary) {
    rlang::abort(
      paste0(
        "shared_strains must contain either 'transmission_score' or ",
        "'is_shared' column."
      ),
      class = "strainflow_column_missing"
    )
  }

  # --- order rows by transmission frequency (most shared at top) ---
  if (has_binary) {
    genome_freq <- shared_strains |>
      dplyr::group_by(.data$genome) |>
      dplyr::summarise(
        freq = sum(.data$is_shared, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$freq)
  } else {
    genome_freq <- shared_strains |>
      dplyr::group_by(.data$genome) |>
      dplyr::summarise(
        freq = sum(.data$transmission_score > 0, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$freq)
  }

  shared_strains <- shared_strains |>
    dplyr::mutate(
      genome = factor(.data$genome, levels = genome_freq$genome)
    )

  # --- determine fill variable ---
  if (has_score) {
    fill_var <- "transmission_score"
  } else {
    shared_strains <- shared_strains |>
      dplyr::mutate(is_shared_num = as.numeric(.data$is_shared))
    fill_var <- "is_shared_num"
  }

  # --- optionally order columns by pair metadata ---
  has_delivery <- !is.null(pair_metadata) &&
    "delivery_mode" %in% colnames(pair_metadata)

  if (has_delivery) {
    pair_order <- pair_metadata |>
      dplyr::arrange(.data$delivery_mode, .data$pair_id) |>
      dplyr::pull(.data$pair_id)
    shared_strains <- shared_strains |>
      dplyr::mutate(
        pair_id = factor(.data$pair_id, levels = pair_order)
      )
  }

  # --- build heatmap ---
  p <- ggplot2::ggplot(
    shared_strains,
    ggplot2::aes(x = .data$pair_id, y = .data$genome,
                 fill = .data[[fill_var]])
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3)

  if (has_score) {
    p <- p +
      viridis::scale_fill_viridis(
        name = "Transmission\nScore",
        option = "D",
        na.value = "grey95"
      )
  } else {
    p <- p +
      ggplot2::scale_fill_gradient(
        low = .sf_colors$not_shared,
        high = .sf_colors$shared_strain,
        name = "Shared",
        labels = c("No", "Yes"),
        breaks = c(0, 1),
        na.value = "grey95"
      )
  }

  p <- p +
    ggplot2::labs(
      title = "Strain Sharing Across Mother-Child Pairs",
      x = "Pair",
      y = "Genome (Species)"
    ) +
    theme_strainflow() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 7)
    )

  # --- add delivery mode annotation strip if available ---
  if (has_delivery) {
    # Build annotation data for the strip
    annot_data <- pair_metadata |>
      dplyr::filter(.data$pair_id %in% levels(shared_strains$pair_id)) |>
      dplyr::mutate(
        pair_id = factor(.data$pair_id, levels = levels(shared_strains$pair_id)),
        y_pos = max(as.numeric(shared_strains$genome), na.rm = TRUE) + 1
      )

    p_annot <- ggplot2::ggplot(
      annot_data,
      ggplot2::aes(x = .data$pair_id, y = 1, fill = .data$delivery_mode)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.3) +
      ggplot2::scale_fill_brewer(palette = "Set2", name = "Delivery Mode") +
      ggplot2::labs(x = NULL, y = NULL) +
      theme_strainflow() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )

    p <- patchwork::wrap_plots(p_annot, p, ncol = 1, heights = c(1, 12))
  }

  .save_or_return(p, output_path, width = 12, height = 10)
}


# ==========================================================================
# Transmission Plots
# ==========================================================================

# --------------------------------------------------------------------------
# 7. plot_transmission_distribution
# --------------------------------------------------------------------------

#' Plot transmission score distribution
#'
#' Overlapping density plots of transmission scores for related (red) vs
#' unrelated (grey) comparisons. A vertical dashed line marks the
#' classification threshold. Individual data points are shown as a rug plot.
#'
#' @param transmission_scores A \code{\link[tibble]{tibble}} with columns
#'   \code{transmission_score} and \code{comparison_type}.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_transmission_distribution <- function(transmission_scores,
                                           output_path = NULL) {

  required <- c("transmission_score", "comparison_type")
  missing_cols <- setdiff(required, colnames(transmission_scores))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "transmission_scores is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  threshold <- 0.99999

  p <- ggplot2::ggplot(
    transmission_scores,
    ggplot2::aes(
      x = .data$transmission_score,
      fill = .data$comparison_type,
      color = .data$comparison_type
    )
  ) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.6) +
    ggplot2::geom_rug(
      ggplot2::aes(color = .data$comparison_type),
      alpha = 0.4,
      sides = "b",
      length = ggplot2::unit(0.03, "npc")
    ) +
    ggplot2::geom_vline(
      xintercept = threshold,
      linetype = "dashed", color = "grey30", linewidth = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = c("related"   = .sf_colors$related,
                 "unrelated" = .sf_colors$unrelated),
      name = "Comparison Type",
      labels = c("related" = "Related", "unrelated" = "Unrelated")
    ) +
    ggplot2::scale_color_manual(
      values = c("related"   = .sf_colors$related,
                 "unrelated" = .sf_colors$unrelated),
      name = "Comparison Type",
      labels = c("related" = "Related", "unrelated" = "Unrelated")
    ) +
    ggplot2::labs(
      title = "Transmission Score Distribution",
      subtitle = glue::glue("Dashed line = threshold ({threshold})"),
      x = "Transmission Score (popANI)",
      y = "Density"
    ) +
    theme_strainflow()

  .save_or_return(p, output_path, width = 9, height = 6)
}


# --------------------------------------------------------------------------
# 8. plot_transmission_rates
# --------------------------------------------------------------------------

#' Plot transmission rate forest plot
#'
#' Forest plot showing species-level transmission rates with confidence
#' intervals. Species are ordered by transmission rate (descending). Significant
#' results are highlighted in red. The number of shared vs tested pairs is
#' annotated on the right margin.
#'
#' @param transmission_rates A \code{\link[tibble]{tibble}} with columns
#'   \code{genome}, \code{transmission_rate}, \code{ci_lower}, \code{ci_upper},
#'   and \code{n_pairs}. Optionally \code{n_shared} for the annotation and
#'   \code{significant} (logical) for coloring.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_transmission_rates <- function(transmission_rates, output_path = NULL) {

  required <- c("genome", "transmission_rate", "ci_lower", "ci_upper", "n_pairs")
  missing_cols <- setdiff(required, colnames(transmission_rates))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "transmission_rates is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  # --- derive n_shared if not present ---
  if (!"n_shared" %in% colnames(transmission_rates)) {
    transmission_rates <- transmission_rates |>
      dplyr::mutate(
        n_shared = round(.data$transmission_rate * .data$n_pairs)
      )
  }

  # --- derive significance if not present ---
  if (!"significant" %in% colnames(transmission_rates)) {
    transmission_rates <- transmission_rates |>
      dplyr::mutate(significant = .data$ci_lower > 0)
  }

  # --- order by transmission rate descending ---
  plot_data <- transmission_rates |>
    dplyr::arrange(.data$transmission_rate) |>
    dplyr::mutate(
      genome = factor(.data$genome, levels = .data$genome),
      sig_label = ifelse(.data$significant, "Significant", "Not significant"),
      annotation = glue::glue("{n_shared}/{n_pairs}")
    )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$transmission_rate,
      y = .data$genome,
      color = .data$sig_label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "solid", color = "grey80", linewidth = 0.4
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
      height = 0.25, linewidth = 0.5
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = max(plot_data$ci_upper, na.rm = TRUE) + 0.05,
        label = .data$annotation
      ),
      hjust = 0, size = 3, color = "grey30"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Significant"     = .sf_colors$significant,
        "Not significant" = .sf_colors$nonsignificant
      ),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(),
      expand = ggplot2::expansion(mult = c(0.02, 0.15))
    ) +
    ggplot2::labs(
      title = "Species-Level Transmission Rates",
      subtitle = "Point estimates with 95% confidence intervals (n shared / n tested)",
      x = "Transmission Rate",
      y = NULL
    ) +
    theme_strainflow() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank()
    )

  # Dynamic height based on number of species
  plot_height <- max(5, nrow(plot_data) * 0.35 + 2)

  .save_or_return(p, output_path, width = 10, height = plot_height)
}


# --------------------------------------------------------------------------
# 9. plot_covariate_comparison
# --------------------------------------------------------------------------

#' Plot covariate comparison
#'
#' Grouped bar chart comparing mean number of shared strains across covariate
#' levels (e.g. vaginal vs cesarean delivery, breast vs formula feeding).
#' Error bars show standard error and significance brackets with p-values are
#' added above the bars.
#'
#' @param covariate_results A \code{\link[tibble]{tibble}} with columns
#'   \code{covariate}, \code{level}, \code{n_shared}, \code{n_total},
#'   \code{rate}, and \code{p_value}.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_covariate_comparison <- function(covariate_results, output_path = NULL) {

  required <- c("covariate", "level", "n_shared", "n_total", "rate", "p_value")
  missing_cols <- setdiff(required, colnames(covariate_results))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "covariate_results is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  # --- compute standard error for error bars ---
  plot_data <- covariate_results |>
    dplyr::mutate(
      se = sqrt(.data$rate * (1 - .data$rate) / pmax(.data$n_total, 1)),
      y_upper = .data$rate + .data$se,
      p_label = dplyr::case_when(
        .data$p_value < 0.001 ~ "***",
        .data$p_value < 0.01  ~ "**",
        .data$p_value < 0.05  ~ "*",
        TRUE                  ~ "ns"
      )
    )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$level, y = .data$rate, fill = .data$level)
  ) +
    ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(.data$rate - .data$se, 0),
                   ymax = .data$rate + .data$se),
      width = 0.2, linewidth = 0.5
    ) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::facet_wrap(~ .data$covariate, scales = "free_x") +
    ggplot2::labs(
      title = "Strain Sharing Rates by Covariate",
      subtitle = "Error bars = standard error; significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant",
      x = NULL,
      y = "Strain Sharing Rate"
    ) +
    theme_strainflow()

  # --- add significance brackets ---
  # For each covariate, add a bracket between the two levels
  bracket_data <- plot_data |>
    dplyr::group_by(.data$covariate) |>
    dplyr::summarise(
      y_bracket = max(.data$y_upper, na.rm = TRUE) * 1.08,
      p_label = dplyr::first(.data$p_label),
      levels_list = list(.data$level),
      .groups = "drop"
    ) |>
    dplyr::filter(purrr::map_int(.data$levels_list, length) >= 2L)

  if (nrow(bracket_data) > 0L) {
    for (i in seq_len(nrow(bracket_data))) {
      row <- bracket_data[i, ]
      lvls <- row$levels_list[[1]]
      p <- p +
        ggplot2::annotate(
          "segment",
          x = 1, xend = min(length(lvls), 2),
          y = row$y_bracket, yend = row$y_bracket,
          linewidth = 0.4
        ) +
        ggplot2::annotate(
          "text",
          x = (1 + min(length(lvls), 2)) / 2,
          y = row$y_bracket * 1.03,
          label = row$p_label,
          size = 3.5, fontface = "bold"
        )
    }
  }

  .save_or_return(p, output_path, width = 10, height = 6)
}


# ==========================================================================
# Longitudinal Plots
# ==========================================================================

# --------------------------------------------------------------------------
# 10. plot_strain_timeline
# --------------------------------------------------------------------------

#' Plot strain timeline per pair
#'
#' Visualizes strain detection and sharing over time for each mother-child
#' pair. Each species is a row, timepoints are on the x-axis, and points
#' indicate strain detection in the mother (triangle) or child (circle).
#' Lines connect mother-child points when the same strain is shared at a
#' given timepoint.
#'
#' @param state_matrix A \code{\link[tibble]{tibble}} from
#'   \code{build_strain_state_matrix} with columns including \code{pair_id},
#'   \code{genome}, \code{timepoint_days}, \code{role}, and \code{detected}
#'   (logical).
#' @param shared_strains Optional \code{\link[tibble]{tibble}} with shared
#'   strain calls containing \code{pair_id}, \code{genome},
#'   \code{timepoint_days}, and \code{is_shared} (logical). Used to draw
#'   connecting lines.
#' @param pair_ids Optional character vector of pair IDs to plot. If
#'   \code{NULL}, all pairs are plotted.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_strain_timeline <- function(state_matrix,
                                 shared_strains = NULL,
                                 pair_ids = NULL,
                                 output_path = NULL) {

  required <- c("pair_id", "genome", "timepoint_days", "role", "detected")
  missing_cols <- setdiff(required, colnames(state_matrix))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "state_matrix is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  # --- filter to requested pairs ---
  plot_data <- state_matrix
  if (!is.null(pair_ids)) {
    plot_data <- plot_data |>
      dplyr::filter(.data$pair_id %in% pair_ids)
  }

  if (nrow(plot_data) == 0L) {
    rlang::warn("No data to plot after filtering by pair_ids.")
    return(ggplot2::ggplot() + theme_strainflow())
  }

  # --- build sharing connection data ---
  connection_data <- NULL
  if (!is.null(shared_strains) &&
      all(c("pair_id", "genome", "timepoint_days", "is_shared") %in%
          colnames(shared_strains))) {

    shared_at_tp <- shared_strains |>
      dplyr::filter(.data$is_shared)

    if (!is.null(pair_ids)) {
      shared_at_tp <- shared_at_tp |>
        dplyr::filter(.data$pair_id %in% pair_ids)
    }

    if (nrow(shared_at_tp) > 0L) {
      # For each shared event, get both mother and child detection points
      mother_pts <- plot_data |>
        dplyr::filter(.data$role == "mother", .data$detected) |>
        dplyr::select("pair_id", "genome", "timepoint_days")

      child_pts <- plot_data |>
        dplyr::filter(.data$role == "child", .data$detected) |>
        dplyr::select("pair_id", "genome", "timepoint_days")

      connection_data <- shared_at_tp |>
        dplyr::select("pair_id", "genome", "timepoint_days") |>
        dplyr::inner_join(mother_pts, by = c("pair_id", "genome", "timepoint_days")) |>
        dplyr::inner_join(child_pts, by = c("pair_id", "genome", "timepoint_days"))
    }
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$timepoint_days,
      y = .data$genome,
      color = .data$genome
    )
  ) +
    # --- detected points ---
    ggplot2::geom_point(
      data = plot_data |> dplyr::filter(.data$detected),
      ggplot2::aes(shape = .data$role),
      size = 3,
      alpha = 0.8
    ) +
    # --- not detected (open circles) ---
    ggplot2::geom_point(
      data = plot_data |> dplyr::filter(!.data$detected),
      ggplot2::aes(shape = .data$role),
      size = 2,
      alpha = 0.3,
      fill = NA
    ) +
    ggplot2::scale_shape_manual(
      values = c("mother" = 17, "child" = 16),
      name = "Role",
      labels = c("mother" = "Mother", "child" = "Child")
    )

  # --- add sharing connection lines ---
  if (!is.null(connection_data) && nrow(connection_data) > 0L) {
    # Draw vertical line segments at shared timepoints between mother and child
    for (i in seq_len(nrow(connection_data))) {
      row <- connection_data[i, ]
      p <- p +
        ggplot2::annotate(
          "segment",
          x = row$timepoint_days, xend = row$timepoint_days,
          y = as.numeric(factor(row$genome, levels = levels(factor(plot_data$genome)))) - 0.15,
          yend = as.numeric(factor(row$genome, levels = levels(factor(plot_data$genome)))) + 0.15,
          color = "grey30",
          linewidth = 0.8,
          alpha = 0.6
        )
    }
  }

  p <- p +
    ggplot2::facet_wrap(~ .data$pair_id, scales = "free_x") +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
    ggplot2::guides(color = "none") +
    ggplot2::labs(
      title = "Strain Detection Timeline",
      subtitle = "Filled shapes = detected; faded = not detected. Triangle = mother, circle = child.",
      x = "Days",
      y = "Genome (Species)"
    ) +
    theme_strainflow() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7)
    )

  n_pairs <- length(unique(plot_data$pair_id))
  plot_height <- max(6, length(unique(plot_data$genome)) * 0.4 + 2)
  plot_width <- max(10, n_pairs * 5)

  .save_or_return(p, output_path, width = plot_width, height = plot_height)
}


# --------------------------------------------------------------------------
# 11. plot_persistence_km
# --------------------------------------------------------------------------

#' Plot persistence Kaplan-Meier curves
#'
#' Plots Kaplan-Meier survival curves for strain persistence after
#' transmission. Includes confidence bands and an at-risk table below the
#' main curve. If a grouping variable is provided, separate curves are
#' drawn for each group.
#'
#' @param persistence_model A list from \code{fit_persistence_model}
#'   containing a \code{km_fit} element (a \code{\link[survival]{survfit}}
#'   object) and optionally \code{surv_data} (the data used to fit the model).
#' @param group_by Optional character scalar. Name of a grouping variable
#'   (e.g. \code{"delivery_mode"}) to stratify the KM curves.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_persistence_km <- function(persistence_model,
                                group_by = NULL,
                                output_path = NULL) {

  if (!is.list(persistence_model) ||
      !"km_fit" %in% names(persistence_model)) {
    rlang::abort(
      paste0(
        "persistence_model must be a list containing a 'km_fit' element ",
        "(a survfit object)."
      ),
      class = "strainflow_invalid_input"
    )
  }

  km_fit <- persistence_model$km_fit

  # --- extract survival data from the survfit object ---
  surv_summary <- summary(km_fit)

  if (is.null(group_by) || is.null(surv_summary$strata)) {
    # Ungrouped KM curve
    surv_df <- tibble::tibble(
      time    = surv_summary$time,
      surv    = surv_summary$surv,
      lower   = surv_summary$lower,
      upper   = surv_summary$upper,
      n.risk  = surv_summary$n.risk,
      n.event = surv_summary$n.event,
      group   = "All"
    )
  } else {
    # Grouped: parse strata names
    strata_names <- names(surv_summary$strata)
    strata_lengths <- surv_summary$strata
    group_vec <- rep(strata_names, times = strata_lengths)
    # Clean strata names: "group_by=value" -> "value"
    group_vec <- stringr::str_replace(group_vec, "^.*=", "")

    surv_df <- tibble::tibble(
      time    = surv_summary$time,
      surv    = surv_summary$surv,
      lower   = surv_summary$lower,
      upper   = surv_summary$upper,
      n.risk  = surv_summary$n.risk,
      n.event = surv_summary$n.event,
      group   = group_vec
    )
  }

  # --- add time 0 row for each group ---
  time_zero <- surv_df |>
    dplyr::distinct(.data$group) |>
    dplyr::mutate(
      time = 0, surv = 1.0, lower = 1.0, upper = 1.0,
      n.risk = NA_real_, n.event = 0
    )
  surv_df <- dplyr::bind_rows(time_zero, surv_df) |>
    dplyr::arrange(.data$group, .data$time)

  # --- main KM plot ---
  p_main <- ggplot2::ggplot(surv_df, ggplot2::aes(x = .data$time, y = .data$surv))

  if (length(unique(surv_df$group)) > 1) {
    p_main <- p_main +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                     fill = .data$group),
        alpha = 0.2
      ) +
      ggplot2::geom_step(
        ggplot2::aes(color = .data$group),
        linewidth = 0.8
      ) +
      ggplot2::scale_color_brewer(palette = "Set1", name = group_by) +
      ggplot2::scale_fill_brewer(palette = "Set1", name = group_by)
  } else {
    p_main <- p_main +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        fill = .sf_colors$reads_after,
        alpha = 0.2
      ) +
      ggplot2::geom_step(
        color = .sf_colors$reads_after,
        linewidth = 0.8
      )
  }

  p_main <- p_main +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(
      title = "Strain Persistence After Transmission",
      subtitle = "Kaplan-Meier survival estimate with 95% confidence bands",
      x = "Days Since Acquisition",
      y = "Probability of Persistence"
    ) +
    theme_strainflow()

  # --- build at-risk table ---
  # Select time points for the table
  time_breaks <- scales::breaks_pretty(n = 6)(range(surv_df$time, na.rm = TRUE))
  time_breaks <- time_breaks[time_breaks >= 0]

  risk_table <- purrr::map_dfr(unique(surv_df$group), function(g) {
    g_data <- surv_df |> dplyr::filter(.data$group == g, !is.na(.data$n.risk))
    purrr::map_dfr(time_breaks, function(t) {
      # Find the closest time point at or before t
      available <- g_data |> dplyr::filter(.data$time <= t)
      if (nrow(available) == 0L) {
        n_at_risk <- NA_real_
      } else {
        n_at_risk <- available |>
          dplyr::filter(.data$time == max(.data$time)) |>
          dplyr::pull(.data$n.risk) |>
          utils::head(1)
      }
      tibble::tibble(group = g, time = t, n_at_risk = n_at_risk)
    })
  })

  p_risk <- ggplot2::ggplot(
    risk_table,
    ggplot2::aes(x = .data$time, y = .data$group,
                 label = as.character(.data$n_at_risk))
  ) +
    ggplot2::geom_text(size = 3, na.rm = TRUE) +
    ggplot2::scale_x_continuous(
      limits = ggplot2::layer_scales(p_main)$x$range$range,
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(x = NULL, y = "At Risk") +
    theme_strainflow() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 9),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(0, 5.5, 5.5, 5.5)
    )

  # --- combine with patchwork ---
  p <- patchwork::wrap_plots(
    p_main + ggplot2::theme(plot.margin = ggplot2::margin(5.5, 5.5, 0, 5.5)),
    p_risk,
    ncol = 1,
    heights = c(4, 1)
  )

  .save_or_return(p, output_path, width = 9, height = 7)
}


# --------------------------------------------------------------------------
# 12. plot_divergence_trajectories
# --------------------------------------------------------------------------

#' Plot post-transmission divergence trajectories
#'
#' Shows how strains diverge after transmission by plotting a chosen metric
#' (popANI or APSS) against time since acquisition. Individual pair-by-species
#' trajectories are shown as thin grey lines with a bold loess trend overlay.
#' If multiple species are present, facets are used.
#'
#' @param divergence A \code{\link[tibble]{tibble}} from
#'   \code{track_divergence} with columns \code{pair_id}, \code{genome},
#'   \code{days_since_acquisition}, and metric columns (\code{popANI} and/or
#'   \code{apss}).
#' @param metric Character scalar. Which metric to plot on the y-axis:
#'   \code{"popANI"} (default) or \code{"apss"}.
#' @param output_path Optional file path to save the plot.
#' @return A ggplot object or file path.
#' @export
plot_divergence_trajectories <- function(divergence,
                                         metric = "popANI",
                                         output_path = NULL) {

  required <- c("pair_id", "genome", "days_since_acquisition")
  missing_cols <- setdiff(required, colnames(divergence))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      glue::glue(
        "divergence is missing required columns: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  metric <- match.arg(metric, choices = c("popANI", "apss"))

  if (!metric %in% colnames(divergence)) {
    rlang::abort(
      glue::glue(
        "Column '{metric}' not found in divergence data. ",
        "Available columns: {paste(colnames(divergence), collapse = ', ')}"
      ),
      class = "strainflow_column_missing"
    )
  }

  # --- create a trajectory identifier ---
  plot_data <- divergence |>
    dplyr::mutate(
      trajectory_id = paste(.data$pair_id, .data$genome, sep = "__")
    )

  n_genomes <- length(unique(plot_data$genome))

  # --- choose y-axis label ---
  y_label <- if (metric == "popANI") {
    "popANI (Population Average Nucleotide Identity)"
  } else {
    "APSS (Average Pairwise Synteny Score)"
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$days_since_acquisition, y = .data[[metric]])
  ) +
    # --- individual trajectories ---
    ggplot2::geom_line(
      ggplot2::aes(group = .data$trajectory_id),
      color = "grey70",
      linewidth = 0.3,
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(group = .data$trajectory_id),
      color = "grey70",
      size = 0.8,
      alpha = 0.4
    )

  # --- add loess smooth ---
  # Use tryCatch in case there are too few data points for loess
  p <- p +
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      color = .sf_colors$reads_after,
      fill = .sf_colors$reads_after,
      linewidth = 1.2,
      alpha = 0.2,
      se = TRUE
    )

  # --- facet by genome if multiple species ---
  if (n_genomes > 1 && n_genomes <= 20) {
    p <- p +
      ggplot2::facet_wrap(~ .data$genome, scales = "free_y")
  } else if (n_genomes > 20) {
    # Too many genomes to facet; color by genome instead
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$days_since_acquisition, y = .data[[metric]])
    ) +
      ggplot2::geom_line(
        ggplot2::aes(group = .data$trajectory_id, color = .data$genome),
        linewidth = 0.3,
        alpha = 0.4
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        formula = y ~ x,
        color = "black",
        linewidth = 1.2,
        alpha = 0.2,
        se = TRUE
      ) +
      ggplot2::guides(color = "none")
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 6),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(
      title = glue::glue("Post-Transmission Divergence: {metric}"),
      subtitle = "Individual trajectories (grey) with loess trend (blue)",
      x = "Days Since Acquisition",
      y = y_label
    ) +
    theme_strainflow()

  plot_width <- if (n_genomes > 1 && n_genomes <= 20) {
    max(10, ceiling(sqrt(n_genomes)) * 4)
  } else {
    10
  }
  plot_height <- if (n_genomes > 1 && n_genomes <= 20) {
    max(7, ceiling(n_genomes / ceiling(sqrt(n_genomes))) * 3.5)
  } else {
    7
  }

  .save_or_return(p, output_path, width = plot_width, height = plot_height)
}
