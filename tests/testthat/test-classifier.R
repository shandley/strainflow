test_that("classify_evolutionary_mode assigns correct categories", {
  # Mock inStrain comparisons with known patterns
  instrain <- tibble::tibble(
    genome = rep(c("species_A", "species_B", "species_C"), each = 10),
    name1 = paste0("sample_", rep(1:10, 3)),
    name2 = paste0("sample_", rep(11:20, 3)),
    popANI = c(
      # species_A: high SNV variance (hypermutator)
      runif(10, 0.9990, 0.99999),
      # species_B: low SNV variance (hyper_recombinator)
      runif(10, 0.99998, 0.999999),
      # species_C: moderate (balanced)
      runif(10, 0.9995, 0.99999)
    ),
    compared_bases_count = rep(100000, 30),
    percent_genome_compared = rep(0.8, 30)
  )

  syntracker <- tibble::tibble(
    genome = rep(c("species_A", "species_B", "species_C"), each = 10),
    sample_id_1 = paste0("sample_", rep(1:10, 3)),
    sample_id_2 = paste0("sample_", rep(11:20, 3)),
    apss = c(
      # species_A: low synteny variance (conserved structure)
      runif(10, 0.95, 0.99),
      # species_B: high synteny variance (hyper_recombinator)
      runif(10, 0.50, 0.95),
      # species_C: moderate
      runif(10, 0.80, 0.95)
    ),
    n_regions_compared = rep(50, 30)
  )

  modes <- classify_evolutionary_mode(instrain, syntracker)
  expect_s3_class(modes, "tbl_df")
  expect_true(all(c("genome", "evolutionary_mode", "popani_weight", "apss_weight") %in% names(modes)))
  expect_true(all(modes$evolutionary_mode %in% c("hypermutator", "hyper_recombinator", "balanced")))
  # Weights should sum to 1
  expect_true(all(abs(modes$popani_weight + modes$apss_weight - 1.0) < 1e-10))
})

test_that("compute_transmission_score handles missing APSS", {
  instrain <- tibble::tibble(
    genome = c("sp_A", "sp_B"),
    name1 = c("s1", "s3"),
    name2 = c("s2", "s4"),
    popANI = c(0.999995, 0.99990),
    compared_bases_count = c(100000, 80000),
    percent_genome_compared = c(0.75, 0.60)
  )

  # Only sp_A has SynTracker data
  syntracker <- tibble::tibble(
    genome = "sp_A",
    sample_id_1 = "s1",
    sample_id_2 = "s2",
    apss = 0.98,
    n_regions_compared = 50
  )

  weights <- tibble::tibble(
    genome = c("sp_A", "sp_B"),
    popani_weight = c(0.5, 0.5),
    apss_weight = c(0.5, 0.5)
  )

  scores <- compute_transmission_score(instrain, syntracker, weights)
  expect_s3_class(scores, "tbl_df")
  expect_equal(nrow(scores), 2)
  # sp_B should still have a score (popANI only)
  expect_true(all(!is.na(scores$transmission_score)))
})

test_that("call_shared_strains applies thresholds correctly", {
  scores <- tibble::tibble(
    genome = c("sp_A", "sp_B", "sp_C"),
    sample_id_1 = rep("s1", 3),
    sample_id_2 = rep("s2", 3),
    popANI = c(0.999999, 0.99990, 0.999995),
    apss = c(0.98, 0.80, NA),
    popani_weight = c(0.5, 0.5, 1.0),
    apss_weight = c(0.5, 0.5, 0.0),
    transmission_score = c(0.999999, 0.99990, 0.999995),
    compared_bases_count = c(100000, 80000, 60000),
    n_regions_compared = c(50, 30, NA)
  )

  shared <- call_shared_strains(scores, threshold = 0.99999)
  expect_true("is_shared" %in% names(shared))
  expect_true("confidence" %in% names(shared))

  # sp_A should be shared (score above threshold)
  expect_true(shared$is_shared[shared$genome == "sp_A"])
  # sp_B should NOT be shared (score below threshold)
  expect_false(shared$is_shared[shared$genome == "sp_B"])
})
