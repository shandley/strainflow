test_that("load_config reads valid YAML", {
  # Create temp config
  tmp <- tempfile(fileext = ".yaml")
  yaml::write_yaml(list(
    project = list(name = "test", output_dir = "output", log_dir = "output/logs"),
    metadata = list(sample_sheet = "samples.csv", pair_sheet = "pairs.csv"),
    tools = list(fastp = NULL, bowtie2 = NULL, samtools = NULL, megahit = NULL, blastn = NULL, conda_env = "test_env"),
    reference = list(human_genome_index = "hg38", species_genome_fasta = "ref.fa",
                     species_genome_index = "ref_bt2", stb_file = "ref.stb"),
    qc = list(threads = 4, min_length = 50, qualified_quality_phred = 20,
              cut_front = TRUE, cut_tail = TRUE, cut_window_size = 4, cut_mean_quality = 20),
    host_removal = list(threads = 8),
    taxonomic_profiling = list(threads = 4, analysis_type = "rel_ab_w_read_stats"),
    read_mapping = list(threads = 8, max_insert_size = 1000),
    assembly = list(threads = 8, min_contig_len = 1000),
    instrain_profile = list(processes = 6, min_read_ani = 0.95, min_mapq = -1,
                            min_cov = 5, min_freq = 0.05, min_genome_coverage = 0.5),
    instrain_compare = list(processes = 6, min_cov = 5, min_freq = 0.05),
    syntracker = list(region_length = 1000, region_spacing = 4000, min_identity = 0.97,
                      min_coverage = 0.70, n_regions = 100, flank_length = 2000),
    transmission = list(popANI_threshold = 0.99999, conANI_threshold = 0.99999,
                        apss_threshold = 0.95, min_compared_bases = 50000,
                        min_compared_regions = 20, min_breadth = 0.5,
                        persistence_window_days = 30, n_negative_controls = 5),
    visualization = list(figure_format = "pdf", figure_dpi = 300,
                         figure_width = 10, figure_height = 8),
    execution = list(workers = 4)
  ), tmp)

  config <- load_config(tmp)
  expect_type(config, "list")
  expect_equal(config$project$name, "test")
  expect_equal(config$qc$threads, 4)
  expect_equal(config$transmission$popANI_threshold, 0.99999)
  unlink(tmp)
})

test_that("load_config errors on missing file", {
  expect_error(load_config("nonexistent.yaml"))
})
