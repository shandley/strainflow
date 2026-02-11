test_that("estimate_strain_prevalence computes background rates", {
  shared <- tibble::tibble(
    genome = rep(c("sp_A", "sp_B"), each = 10),
    sample_id_1 = paste0("s", 1:20),
    sample_id_2 = paste0("s", 21:40),
    pair_id = rep(paste0("ctrl_", 1:10), 2),
    comparison_type = rep("negative_control", 20),
    is_shared = c(rep(TRUE, 2), rep(FALSE, 8), rep(TRUE, 1), rep(FALSE, 9)),
    transmission_score = c(rep(0.999999, 2), rep(0.999, 8), rep(0.999999, 1), rep(0.999, 9))
  )

  sample_meta <- tibble::tibble(
    sample_id = c(paste0("s", 1:40)),
    subject_id = c(paste0("subj", 1:40)),
    role = rep(c("mother", "child"), 20),
    timepoint_days = rep(0, 40),
    pair_id = c(rep(paste0("ctrl_", 1:10), 2), rep(paste0("ctrl_", 1:10), 2))
  )

  pair_meta <- tibble::tibble(
    pair_id = paste0("P", 1:3),
    mother_subject_id = paste0("M", 1:3),
    child_subject_id = paste0("C", 1:3),
    delivery_mode = c("vaginal", "cesarean", "vaginal")
  )

  prev <- estimate_strain_prevalence(shared, sample_meta, pair_meta)
  expect_s3_class(prev, "tbl_df")
  expect_true(all(c("genome", "background_sharing_rate") %in% names(prev)))
  # sp_A: 2/10 = 0.2, sp_B: 1/10 = 0.1
  expect_equal(prev$background_sharing_rate[prev$genome == "sp_A"], 0.2)
  expect_equal(prev$background_sharing_rate[prev$genome == "sp_B"], 0.1)
})

test_that("classify_transmission_event assigns correct classes", {
  odds <- tibble::tibble(
    genome = c("sp_A", "sp_B", "sp_C"),
    pair_id = rep("P001", 3),
    sample_id_mother = rep("M001_T0", 3),
    sample_id_child = rep("C001_T0", 3),
    timepoint_days = rep(0, 3),
    observed_sharing = c(TRUE, TRUE, FALSE),
    background_rate = c(0.01, 0.30, 0.01),
    transmission_odds = c(150, 1.5, 0.1),
    log_odds = log10(c(150, 1.5, 0.1))
  )

  classified <- classify_transmission_event(odds)
  expect_true(all(c("transmission_class", "evidence_summary") %in% names(classified)))
  expect_equal(classified$transmission_class[classified$genome == "sp_A"], "high_confidence")
  expect_equal(classified$transmission_class[classified$genome == "sp_B"], "plausible")
  expect_equal(classified$transmission_class[classified$genome == "sp_C"], "unlikely")
})
