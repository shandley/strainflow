# _targets.R — StrainFlow Pipeline Definition
# Run with: targets::tar_make()
# Visualize with: targets::tar_visnetwork()

library(targets)
library(tarchetypes)

# Source all R functions
tar_source("R")

# Optional: parallel execution with crew
# library(crew)
# tar_option_set(
#   controller = crew_controller_local(workers = 4, seconds_idle = 60)
# )

list(
  # ============================================================
  # CONFIGURATION AND METADATA
  # ============================================================

  tar_target(config_file, "config.yaml", format = "file"),
  tar_target(config, load_config(config_file)),

  tar_target(
    sample_sheet_file,
    config$metadata$sample_sheet,
    format = "file"
  ),
  tar_target(sample_metadata, load_sample_metadata(sample_sheet_file)),

  tar_target(
    pair_sheet_file,
    config$metadata$pair_sheet,
    format = "file"
  ),
  tar_target(pair_metadata, load_pair_metadata(pair_sheet_file)),

  tar_target(output_dirs, ensure_output_dirs(config)),
  tar_target(tools_ok, validate_tool_paths(config)),

  # ============================================================
  # STEP 1: QC WITH FASTP
  # ============================================================

  tar_target(
    qc_results,
    run_fastp(
      sample_id  = sample_metadata$sample_id,
      fastq_r1   = sample_metadata$fastq_r1,
      fastq_r2   = sample_metadata$fastq_r2,
      output_dir = fs::path(config$project$output_dir, "01_qc"),
      config     = config
    ),
    pattern = map(sample_metadata),
    iteration = "list"
  ),

  tar_target(
    qc_stats,
    aggregate_qc_stats(dplyr::bind_rows(qc_results))
  ),

  # ============================================================
  # STEP 2: HOST READ REMOVAL
  # ============================================================

  tar_target(
    host_results,
    run_host_removal(
      sample_id  = qc_results[[1]]$sample_id,
      trimmed_r1 = qc_results[[1]]$trimmed_r1,
      trimmed_r2 = qc_results[[1]]$trimmed_r2,
      output_dir = fs::path(config$project$output_dir, "02_host_removed"),
      config     = config
    ),
    pattern = map(qc_results),
    iteration = "list"
  ),

  # ============================================================
  # STEP 3: TAXONOMIC PROFILING
  # ============================================================

  tar_target(
    metaphlan_results,
    run_metaphlan4(
      sample_id  = host_results[[1]]$sample_id,
      nonhost_r1 = host_results[[1]]$nonhost_r1,
      nonhost_r2 = host_results[[1]]$nonhost_r2,
      output_dir = fs::path(config$project$output_dir, "03_taxonomic_profiles"),
      config     = config
    ),
    pattern = map(host_results),
    iteration = "list"
  ),

  tar_target(
    parsed_profiles,
    purrr::map2_dfr(
      purrr::map_chr(metaphlan_results, "profile_path"),
      purrr::map_chr(metaphlan_results, "sample_id"),
      parse_metaphlan_profile
    )
  ),

  tar_target(merged_profiles, merge_metaphlan_profiles(parsed_profiles)),

  # ============================================================
  # STEP 4: READ MAPPING TO REFERENCE GENOMES
  # ============================================================

  tar_target(
    mapping_results,
    run_read_mapping(
      sample_id  = host_results[[1]]$sample_id,
      nonhost_r1 = host_results[[1]]$nonhost_r1,
      nonhost_r2 = host_results[[1]]$nonhost_r2,
      output_dir = fs::path(config$project$output_dir, "04_mapped_reads"),
      config     = config
    ),
    pattern = map(host_results),
    iteration = "list"
  ),

  # ============================================================
  # STEP 5: METAGENOME ASSEMBLY (for SynTracker)
  # ============================================================

  tar_target(
    assembly_results,
    run_megahit(
      sample_id  = host_results[[1]]$sample_id,
      nonhost_r1 = host_results[[1]]$nonhost_r1,
      nonhost_r2 = host_results[[1]]$nonhost_r2,
      output_dir = fs::path(config$project$output_dir, "05_assembly"),
      config     = config
    ),
    pattern = map(host_results),
    iteration = "list"
  ),

  # ============================================================
  # STEP 6: INSTRAIN (SNV SIGNAL)
  # ============================================================

  tar_target(
    instrain_profiles,
    run_instrain_profile(
      sample_id  = mapping_results[[1]]$sample_id,
      bam_path   = mapping_results[[1]]$bam_path,
      output_dir = fs::path(config$project$output_dir, "06_instrain"),
      config     = config
    ),
    pattern = map(mapping_results),
    iteration = "list"
  ),

  tar_target(
    genome_info_all,
    purrr::map2_dfr(
      purrr::map_chr(instrain_profiles, "genome_info_path"),
      purrr::map_chr(instrain_profiles, "sample_id"),
      parse_instrain_genome_info
    )
  ),

  # Group profiles by pair for inStrain compare
  tar_target(
    pair_profile_groups,
    build_pair_profile_groups(
      dplyr::bind_rows(instrain_profiles),
      sample_metadata,
      pair_metadata
    )
  ),

  tar_target(
    instrain_compare_results,
    run_instrain_compare(
      comparison_group_id = pair_profile_groups$pair_id,
      is_profile_dirs     = pair_profile_groups$profile_dirs,
      output_dir = fs::path(config$project$output_dir, "06_instrain"),
      config     = config
    ),
    pattern = map(pair_profile_groups),
    iteration = "list"
  ),

  tar_target(
    all_instrain_comparisons,
    purrr::map_dfr(
      purrr::map_chr(instrain_compare_results, "comparisons_table_path"),
      parse_instrain_comparisons
    )
  ),

  # Negative controls
  tar_target(
    negative_manifest,
    generate_negative_controls(sample_metadata, pair_metadata,
                               n_controls = config$transmission$n_negative_controls)
  ),

  # ============================================================
  # STEP 7: SYNTRACKER (SYNTENY SIGNAL)
  # ============================================================

  tar_target(
    syntracker_results,
    run_syntracker_analysis(
      assembly_results = dplyr::bind_rows(assembly_results),
      reference_fasta  = config$reference$species_genome_fasta,
      sample_metadata  = sample_metadata,
      pair_metadata    = pair_metadata,
      config           = config
    )
  ),

  # ============================================================
  # STEP 8: MULTI-SIGNAL TRANSMISSION CLASSIFIER
  # ============================================================

  tar_target(
    evolutionary_modes,
    classify_evolutionary_mode(
      all_instrain_comparisons,
      syntracker_results
    )
  ),

  tar_target(
    species_weights,
    train_species_weights(
      all_instrain_comparisons,
      syntracker_results,
      sample_metadata,
      evolutionary_modes
    )
  ),

  tar_target(
    transmission_scores,
    compute_transmission_score(
      all_instrain_comparisons,
      syntracker_results,
      species_weights
    )
  ),

  tar_target(
    shared_strains,
    call_shared_strains(
      transmission_scores,
      threshold    = config$transmission$popANI_threshold,
      min_compared = config$transmission$min_compared_bases,
      min_regions  = config$transmission$min_compared_regions
    )
  ),

  # ============================================================
  # STEP 9: TRANSMISSION VS ENVIRONMENTAL ACQUISITION
  # ============================================================

  tar_target(
    strain_prevalence,
    estimate_strain_prevalence(shared_strains, sample_metadata, pair_metadata)
  ),

  tar_target(
    transmission_odds,
    compute_transmission_odds(shared_strains, strain_prevalence, sample_metadata)
  ),

  tar_target(
    transmission_model,
    fit_transmission_model(shared_strains, sample_metadata, pair_metadata, strain_prevalence)
  ),

  tar_target(
    classified_events,
    classify_transmission_event(transmission_odds, transmission_model)
  ),

  # ============================================================
  # STEP 10: LONGITUDINAL STRAIN DYNAMICS
  # ============================================================

  tar_target(
    strain_state_matrix,
    build_strain_state_matrix(shared_strains, sample_metadata)
  ),

  tar_target(
    acquisition_events,
    detect_acquisition_events(strain_state_matrix, sample_metadata)
  ),

  tar_target(
    persistence_model,
    fit_persistence_model(strain_state_matrix, acquisition_events, pair_metadata)
  ),

  tar_target(
    divergence_data,
    track_divergence(
      all_instrain_comparisons,
      syntracker_results,
      acquisition_events,
      sample_metadata
    )
  ),

  tar_target(
    replacement_events,
    detect_strain_replacement(strain_state_matrix, all_instrain_comparisons, sample_metadata)
  ),

  # ============================================================
  # STEP 11: VISUALIZATION
  # ============================================================

  tar_target(
    fig_qc,
    plot_qc_summary(
      qc_stats, sample_metadata,
      output_path = fs::path(config$project$output_dir, "09_figures", "01_qc_summary.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_host,
    plot_host_removal(
      dplyr::bind_rows(host_results),
      output_path = fs::path(config$project$output_dir, "09_figures", "02_host_removal.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_taxonomy,
    plot_taxonomic_composition(
      merged_profiles$long_format, sample_metadata,
      output_path = fs::path(config$project$output_dir, "09_figures", "03_taxonomy.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_popani_apss,
    plot_popani_vs_apss(
      transmission_scores,
      output_path = fs::path(config$project$output_dir, "09_figures", "04_popani_vs_apss.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_evo_modes,
    plot_evolutionary_modes(
      evolutionary_modes,
      output_path = fs::path(config$project$output_dir, "09_figures", "05_evolutionary_modes.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_heatmap,
    plot_strain_sharing_heatmap(
      shared_strains, pair_metadata,
      output_path = fs::path(config$project$output_dir, "09_figures", "06_sharing_heatmap.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_transmission_dist,
    plot_transmission_distribution(
      transmission_scores,
      output_path = fs::path(config$project$output_dir, "09_figures", "07_transmission_dist.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_rates,
    plot_transmission_rates(
      classified_events,
      output_path = fs::path(config$project$output_dir, "09_figures", "08_transmission_rates.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_covariates,
    plot_covariate_comparison(
      classified_events,
      output_path = fs::path(config$project$output_dir, "09_figures", "09_covariates.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_timeline,
    plot_strain_timeline(
      strain_state_matrix, shared_strains,
      output_path = fs::path(config$project$output_dir, "09_figures", "10_strain_timeline.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_km,
    plot_persistence_km(
      persistence_model,
      output_path = fs::path(config$project$output_dir, "09_figures", "11_persistence_km.pdf")
    ),
    format = "file"
  ),

  tar_target(
    fig_divergence,
    plot_divergence_trajectories(
      divergence_data,
      output_path = fs::path(config$project$output_dir, "09_figures", "12_divergence.pdf")
    ),
    format = "file"
  )
)
