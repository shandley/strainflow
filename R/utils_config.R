# ============================================================================
# utils_config.R - Pipeline configuration loading, validation, and setup
# ============================================================================

#' Required top-level keys in the pipeline configuration
#' @noRd
.config_required_keys <- c(

"project", "metadata", "tools", "reference", "qc", "host_removal",
  "taxonomic_profiling", "read_mapping", "instrain_profile",
  "instrain_compare", "syntracker", "transmission", "visualization",
  "execution"
)

#' Standard output subdirectories for the pipeline
#' @noRd
.output_subdirs <- c(
  "01_qc",
  "02_host_removed",
  "03_taxonomic_profiles",

  "04_mapped_reads",
  "05_assembly",
  "06_instrain",
  "07_syntracker",
  "08_transmission",
  "09_figures",
  "logs"
)

#' Load and validate pipeline configuration
#'
#' Reads a YAML configuration file, validates that all required top-level keys
#' are present, and fills in sensible defaults for optional parameters.
#'
#' @param config_path Path to config.yaml. Defaults to \code{"config.yaml"} in
#'   the current working directory.
#' @return Named list of validated configuration parameters.
#' @export
load_config <- function(config_path = "config.yaml") {
  # --- check file existence ---
  if (!fs::file_exists(config_path)) {
    rlang::abort(
      glue::glue("Configuration file not found: {config_path}"),
      class = "strainflow_config_missing"
    )
  }

  config <- tryCatch(
    yaml::yaml.load_file(config_path),
    error = function(e) {
      rlang::abort(
        glue::glue(
          "Failed to parse YAML configuration at {config_path}: {conditionMessage(e)}"
        ),
        class = "strainflow_config_parse_error",
        parent = e
      )
    }
  )

  if (!is.list(config)) {
    rlang::abort(
      "Configuration file did not parse to a named list. Check YAML format.",
      class = "strainflow_config_format_error"
    )
  }

  # --- validate required keys ---
  present_keys <- names(config)
  missing_keys <- setdiff(.config_required_keys, present_keys)

  if (length(missing_keys) > 0L) {
    rlang::abort(
      glue::glue(
        "Configuration is missing required top-level keys: ",
        "{paste(missing_keys, collapse = ', ')}"
      ),
      class = "strainflow_config_missing_keys"
    )
  }

  # --- apply defaults ---
  config <- .apply_config_defaults(config)

  # --- resolve output_dir to absolute path ---
  if (!is.null(config$project$output_dir)) {
    config$project$output_dir <- fs::path_abs(config$project$output_dir)
  }

  cli::cli_inform(c(
    "v" = "Configuration loaded from {.file {config_path}}",
    "i" = "Project: {.val {config$project$name %||% 'unnamed'}}"
  ))

  config
}

#' Apply sensible defaults for optional config parameters
#'
#' @param config Raw parsed config list.
#' @return Config list with defaults filled in.
#' @noRd
.apply_config_defaults <- function(config) {
  # -- project defaults --
  config$project$name      <- config$project$name      %||% "strainflow_run"
  config$project$output_dir <- config$project$output_dir %||% "output"

  # -- qc defaults (fastp) --
  config$qc$min_quality   <- config$qc$min_quality   %||% 20L

config$qc$min_length    <- config$qc$min_length    %||% 50L
  config$qc$adapter_trim  <- config$qc$adapter_trim  %||% TRUE
  config$qc$threads       <- config$qc$threads       %||% 4L

  # -- host_removal defaults --
  config$host_removal$threads <- config$host_removal$threads %||% 4L

  # -- read_mapping defaults --
  config$read_mapping$threads   <- config$read_mapping$threads   %||% 8L
  config$read_mapping$min_mapq  <- config$read_mapping$min_mapq  %||% 2L

  # -- instrain_profile defaults --
  config$instrain_profile$min_coverage    <-
    config$instrain_profile$min_coverage    %||% 5
  config$instrain_profile$min_freq        <-
    config$instrain_profile$min_freq        %||% 0.05
  config$instrain_profile$fdr             <-
    config$instrain_profile$fdr             %||% 1e-6
  config$instrain_profile$min_genome_cov  <-
    config$instrain_profile$min_genome_cov  %||% 0.5

  # -- instrain_compare defaults --
  config$instrain_compare$min_coverage  <-
    config$instrain_compare$min_coverage  %||% 5
  config$instrain_compare$min_freq      <-
    config$instrain_compare$min_freq      %||% 0.05

  # -- syntracker defaults --
  config$syntracker$synteny_threshold <-
    config$syntracker$synteny_threshold %||% 0.99

  # -- transmission classifier defaults --
  config$transmission$popani_threshold      <-
    config$transmission$popani_threshold      %||% 0.99999
  config$transmission$synteny_threshold     <-
    config$transmission$synteny_threshold     %||% 0.99
  config$transmission$min_shared_bases      <-
    config$transmission$min_shared_bases      %||% 10000L
  config$transmission$negative_control_fdr  <-
    config$transmission$negative_control_fdr  %||% 0.05

  # -- execution defaults --
  config$execution$workers    <- config$execution$workers    %||% 4L
  config$execution$memory_gb  <- config$execution$memory_gb  %||% 16L

  config
}

#' Validate external tool availability
#'
#' Checks that every external bioinformatics tool required by the pipeline is
#' either found on the system \code{PATH} or at the location specified in the
#' pipeline configuration. Also verifies that the conda environments for
#' inStrain and MetaPhlAn4 exist.
#'
#' @param config Pipeline config list as returned by \code{\link{load_config}}.
#' @return Invisible \code{TRUE} on success. Stops with an informative error if
#'   any tools are missing.
#' @export
validate_tool_paths <- function(config) {
  required_tools <- c("fastp", "bowtie2", "samtools", "megahit")
  missing <- character()

  for (tool in required_tools) {
    resolved <- tryCatch(
      resolve_tool(tool, config),
      error = function(e) NULL
    )
    if (is.null(resolved)) {
      missing <- c(missing, tool)
    }
  }

  if (length(missing) > 0L) {
    rlang::abort(
      c(
        glue::glue(
          "The following required tools could not be found: ",
          "{paste(missing, collapse = ', ')}"
        ),
        "i" = "Ensure each tool is on your PATH or specified under `tools:` in config.yaml."
      ),
      class = "strainflow_tools_missing"
    )
  }

  # -- check conda environments for inStrain and MetaPhlAn4 --
  conda_envs_required <- c(
    instrain = config$tools$instrain_env %||% "instrain",
    metaphlan = config$tools$metaphlan_env %||% "metaphlan4"
  )

  conda_list_result <- tryCatch(
    processx::run("conda", args = c("env", "list"), error_on_status = FALSE),
    error = function(e) {
      rlang::abort(
        c(
          "Could not run `conda env list`.",
          "i" = "Is conda installed and on your PATH?"
        ),
        class = "strainflow_conda_not_found",
        parent = e
      )
    }
  )

  env_lines <- conda_list_result$stdout
  missing_envs <- character()

  for (nm in names(conda_envs_required)) {
    env_name <- conda_envs_required[[nm]]
    # Match env name at start of line (conda env list format)
    pattern <- paste0("(?m)^", env_name, "\\s")
    if (!stringr::str_detect(env_lines, pattern)) {
      missing_envs <- c(missing_envs, glue::glue("{nm} (env: {env_name})"))
    }
  }

  if (length(missing_envs) > 0L) {
    rlang::abort(
      c(
        glue::glue(
          "The following conda environments are missing: ",
          "{paste(missing_envs, collapse = ', ')}"
        ),
        "i" = "Create the environments or update the env names under `tools:` in config.yaml."
      ),
      class = "strainflow_conda_envs_missing"
    )
  }

  cli::cli_inform(c("v" = "All required tools and conda environments validated."))
  invisible(TRUE)
}

#' Resolve tool path from config or PATH
#'
#' Looks up a tool first in the \code{config$tools} list, then falls back to
#' \code{Sys.which()} for PATH resolution. Returns the absolute path to the
#' tool executable.
#'
#' @param tool_name Tool name string (e.g. \code{"fastp"}).
#' @param config Pipeline config list.
#' @return Absolute path to the tool executable.
#' @noRd
resolve_tool <- function(tool_name, config) {
  # 1. Check config$tools for an explicit path
  config_path <- config$tools[[tool_name]]

  if (!is.null(config_path) && nzchar(config_path)) {
    abs_path <- fs::path_abs(config_path)
    if (fs::file_exists(abs_path)) {
      return(as.character(abs_path))
    }
    rlang::abort(
      glue::glue(
        "Tool '{tool_name}' specified in config at '{config_path}' ",
        "but the file does not exist (resolved to: {abs_path})."
      ),
      class = "strainflow_tool_not_found"
    )
  }

  # 2. Fall back to PATH
  sys_path <- Sys.which(tool_name)
  if (nzchar(sys_path)) {
    return(as.character(fs::path_abs(sys_path)))
  }

  rlang::abort(
    glue::glue(
      "Tool '{tool_name}' is not specified in config and was not found on PATH."
    ),
    class = "strainflow_tool_not_found"
  )
}

#' Create output directory structure
#'
#' Creates the full set of numbered output subdirectories plus a \code{logs/}
#' directory under the project output root specified in the config.
#'
#' @param config Pipeline config list as returned by \code{\link{load_config}}.
#' @return Character vector of created directory paths (invisibly).
#' @export
ensure_output_dirs <- function(config) {
  output_root <- config$project$output_dir

  if (is.null(output_root) || !nzchar(output_root)) {
    rlang::abort(
      "config$project$output_dir is not set. Cannot create output directories.",
      class = "strainflow_config_missing_output"
    )
  }

  dirs_to_create <- fs::path(output_root, .output_subdirs)

  purrr::walk(dirs_to_create, function(d) {
    fs::dir_create(d, recurse = TRUE)
  })

  cli::cli_inform(c(
    "v" = "Output directory structure created under {.path {output_root}}",
    "i" = "{length(dirs_to_create)} subdirectories ensured."
  ))

  invisible(as.character(dirs_to_create))
}
