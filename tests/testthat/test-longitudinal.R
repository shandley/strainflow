test_that("build_strain_state_matrix creates correct structure", {
  shared <- tibble::tibble(
    genome = rep("sp_A", 4),
    sample_id_1 = c("M001_T0", "M001_T0", "M001_T30", "M001_T30"),
    sample_id_2 = c("C001_T0", "C001_T30", "C001_T0", "C001_T30"),
    is_shared = c(TRUE, TRUE, TRUE, FALSE),
    transmission_score = c(0.999999, 0.999998, 0.999997, 0.99990)
  )

  sample_meta <- tibble::tibble(
    sample_id = c("M001_T0", "M001_T30", "C001_T0", "C001_T30"),
    subject_id = c("M001", "M001", "C001", "C001"),
    role = c("mother", "mother", "child", "child"),
    timepoint_days = c(0, 30, 0, 30),
    pair_id = rep("P001", 4)
  )

  state <- build_strain_state_matrix(shared, sample_meta)
  expect_s3_class(state, "tbl_df")
  expect_true(all(c("pair_id", "genome", "subject_id", "role",
                     "timepoint_days", "strain_detected") %in% names(state)))
})

test_that("detect_acquisition_events classifies timing correctly", {
  state <- tibble::tibble(
    pair_id = rep("P001", 6),
    genome = rep(c("sp_A", "sp_B"), each = 3),
    subject_id = rep("C001", 6),
    role = rep("child", 6),
    timepoint_days = rep(c(0, 30, 90), 2),
    strain_detected = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE)
  )

  sample_meta <- tibble::tibble(
    sample_id = c("M001_T0", "C001_T0", "C001_T30", "C001_T90"),
    subject_id = c("M001", "C001", "C001", "C001"),
    role = c("mother", "child", "child", "child"),
    timepoint_days = c(0, 0, 30, 90),
    pair_id = rep("P001", 4)
  )

  events <- detect_acquisition_events(state, sample_meta)
  expect_s3_class(events, "tbl_df")

  # sp_A: detected at day 0 → "birth"
  sp_a <- events[events$genome == "sp_A", ]
  expect_equal(sp_a$acquisition_timepoint_days, 0)
  expect_equal(sp_a$acquisition_type, "birth")

  # sp_B: detected at day 90 → "delayed"
  sp_b <- events[events$genome == "sp_B", ]
  expect_equal(sp_b$acquisition_timepoint_days, 90)
  expect_equal(sp_b$acquisition_type, "delayed")
})

test_that("detect_strain_replacement identifies replacement events", {
  state <- tibble::tibble(
    pair_id = rep("P001", 3),
    genome = rep("sp_A", 3),
    subject_id = rep("C001", 3),
    role = rep("child", 3),
    timepoint_days = c(0, 30, 90),
    strain_detected = c(TRUE, TRUE, TRUE)
  )

  # Within-child comparisons show strain change at T90
  instrain <- tibble::tibble(
    sample_id_1 = c("C001_T0", "C001_T0", "C001_T30"),
    sample_id_2 = c("C001_T30", "C001_T90", "C001_T90"),
    genome = rep("sp_A", 3),
    popANI = c(0.999999, 0.99950, 0.99950),
    compared_bases_count = rep(100000, 3),
    percent_genome_compared = rep(0.8, 3)
  )

  sample_meta <- tibble::tibble(
    sample_id = c("C001_T0", "C001_T30", "C001_T90"),
    subject_id = rep("C001", 3),
    role = rep("child", 3),
    timepoint_days = c(0, 30, 90),
    pair_id = rep("P001", 3)
  )

  replacements <- detect_strain_replacement(state, instrain, sample_meta)
  expect_s3_class(replacements, "tbl_df")
  # Should detect replacement between T30 and T90 (popANI dropped)
  expect_true(nrow(replacements) > 0)
})
