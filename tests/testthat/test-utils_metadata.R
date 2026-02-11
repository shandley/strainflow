test_that("load_sample_metadata validates required columns", {
  tmp <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(
    sample_id = c("M001_T0", "C001_T0"),
    subject_id = c("M001", "C001"),
    role = c("mother", "child"),
    timepoint_days = c(0, 0),
    pair_id = c("P001", "P001"),
    fastq_r1 = c("r1.fq.gz", "r1.fq.gz"),
    fastq_r2 = c("r2.fq.gz", "r2.fq.gz")
  ), tmp)

  meta <- load_sample_metadata(tmp)
  expect_s3_class(meta, "tbl_df")
  expect_equal(nrow(meta), 2)
  expect_true(all(c("sample_id", "subject_id", "role", "timepoint_days", "pair_id") %in% names(meta)))
  unlink(tmp)
})

test_that("load_sample_metadata errors on missing columns", {
  tmp <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(
    sample_id = "M001",
    subject_id = "M001"
  ), tmp)

  expect_error(load_sample_metadata(tmp))
  unlink(tmp)
})

test_that("load_pair_metadata works", {
  tmp <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(
    pair_id = "P001",
    mother_subject_id = "M001",
    child_subject_id = "C001",
    delivery_mode = "vaginal"
  ), tmp)

  pairs <- load_pair_metadata(tmp)
  expect_s3_class(pairs, "tbl_df")
  expect_equal(nrow(pairs), 1)
  unlink(tmp)
})

test_that("generate_comparison_manifest creates correct pairs", {
  sample_meta <- tibble::tibble(
    sample_id = c("M001_T0", "M001_T30", "C001_T0", "C001_T30"),
    subject_id = c("M001", "M001", "C001", "C001"),
    role = c("mother", "mother", "child", "child"),
    timepoint_days = c(0, 30, 0, 30),
    pair_id = c("P001", "P001", "P001", "P001"),
    fastq_r1 = rep("r1.fq", 4),
    fastq_r2 = rep("r2.fq", 4)
  )

  pair_meta <- tibble::tibble(
    pair_id = "P001",
    mother_subject_id = "M001",
    child_subject_id = "C001",
    delivery_mode = "vaginal"
  )

  manifest <- generate_comparison_manifest(sample_meta, pair_meta)
  expect_s3_class(manifest, "tbl_df")

  # Should have mother_child, within_mother, within_child comparisons
  expect_true("mother_child" %in% manifest$comparison_type)
  # 2 mother x 2 child = 4 mother_child comparisons
  expect_equal(sum(manifest$comparison_type == "mother_child"), 4)
})

test_that("generate_negative_controls excludes related pairs", {
  sample_meta <- tibble::tibble(
    sample_id = c("M001_T0", "C001_T0", "M002_T0", "C002_T0"),
    subject_id = c("M001", "C001", "M002", "C002"),
    role = c("mother", "child", "mother", "child"),
    timepoint_days = c(0, 0, 0, 0),
    pair_id = c("P001", "P001", "P002", "P002"),
    fastq_r1 = rep("r1.fq", 4),
    fastq_r2 = rep("r2.fq", 4)
  )

  pair_meta <- tibble::tibble(
    pair_id = c("P001", "P002"),
    mother_subject_id = c("M001", "M002"),
    child_subject_id = c("C001", "C002"),
    delivery_mode = c("vaginal", "cesarean")
  )

  controls <- generate_negative_controls(sample_meta, pair_meta, n_controls = 2)
  expect_s3_class(controls, "tbl_df")
  expect_true(all(controls$comparison_type == "negative_control"))
})
