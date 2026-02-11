# ============================================================================
# utils_shell.R - External command execution and conda environment management
# ============================================================================

#' Run an external command with processx
#'
#' Executes an external program, captures standard output and standard error,
#' optionally writes combined output to a log file, and raises an informative
#' error on non-zero exit.
#'
#' @param cmd Path to the executable (resolved via \code{Sys.which} if not
#'   absolute).
#' @param args Character vector of command-line arguments. Defaults to an
#'   empty vector.
#' @param log_file Optional path to a log file. If provided, stdout and stderr
#'   are written to this file.
#' @param timeout_sec Timeout in seconds. Defaults to 7200 (2 hours).
#' @return A named list with elements \code{status} (integer exit code),
#'   \code{stdout} (character), and \code{stderr} (character).
#' @export
run_cmd <- function(cmd,
                    args = character(),
                    log_file = NULL,
                    timeout_sec = 7200) {
  # -- validate cmd exists --
  if (!fs::file_exists(cmd) && !nzchar(Sys.which(cmd))) {
    rlang::abort(
      glue::glue("Command not found: {cmd}"),
      class = "strainflow_cmd_not_found"
    )
  }

  cmd_display <- paste(c(cmd, args), collapse = " ")

  cli::cli_inform(c(
    "i" = "Running: {.code {cmd_display}}"
  ))

  result <- tryCatch(
    processx::run(
      command            = cmd,
      args               = args,
      error_on_status    = FALSE,
      timeout            = timeout_sec,
      spinner            = FALSE,
      echo_cmd           = FALSE,
      stderr_to_stdout   = FALSE
    ),
    error = function(e) {
      rlang::abort(
        c(
          glue::glue("Failed to execute command: {cmd_display}"),
          "x" = conditionMessage(e)
        ),
        class = "strainflow_cmd_exec_error",
        parent = e
      )
    }
  )

  # -- write log file --
  if (!is.null(log_file)) {
    fs::dir_create(fs::path_dir(log_file), recurse = TRUE)

    log_content <- paste0(
      "=== STRAINFLOW COMMAND LOG ===\n",
      "Command: ", cmd_display, "\n",
      "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      "Exit code: ", result$status, "\n",
      "Timeout: ", result$timeout %||% FALSE, "\n",
      "\n=== STDOUT ===\n",
      result$stdout,
      "\n=== STDERR ===\n",
      result$stderr,
      "\n"
    )

    writeLines(log_content, con = log_file)
  }

  # -- check for timeout --
  if (isTRUE(result$timeout)) {
    rlang::abort(
      c(
        glue::glue("Command timed out after {timeout_sec} seconds."),
        "x" = glue::glue("Command: {cmd_display}"),
        "i" = if (!is.null(log_file)) {
          glue::glue("See log: {log_file}")
        } else {
          "Set log_file to capture output for debugging."
        }
      ),
      class = "strainflow_cmd_timeout"
    )
  }

  # -- check for non-zero exit --
  if (result$status != 0L) {
    stderr_lines <- strsplit(result$stderr, "\n", fixed = TRUE)[[1]]
    n_lines <- length(stderr_lines)
    tail_lines <- utils::tail(stderr_lines, 20L)
    tail_display <- paste(tail_lines, collapse = "\n")

    rlang::abort(
      c(
        glue::glue(
          "Command exited with status {result$status}: {cmd_display}"
        ),
        "x" = glue::glue(
          "Last {min(20L, n_lines)} lines of stderr:\n{tail_display}"
        ),
        "i" = if (!is.null(log_file)) {
          glue::glue("Full log: {log_file}")
        } else {
          "Set log_file to capture full output."
        }
      ),
      class = "strainflow_cmd_failed"
    )
  }

  cli::cli_inform(c(
    "v" = "Command completed successfully (exit 0)."
  ))

  list(
    status = result$status,
    stdout = result$stdout,
    stderr = result$stderr
  )
}

#' Run command in conda environment
#'
#' Activates a named conda environment and runs a command inside it. This is
#' the primary mechanism for calling Python-based tools (inStrain, MetaPhlAn4,
#' SynTracker) from R.
#'
#' @param conda_env Conda environment name (e.g. \code{"instrain"}).
#' @param cmd Command to run inside the activated environment.
#' @param args Character vector of arguments. Defaults to an empty vector.
#' @param ... Additional arguments passed to \code{\link{run_cmd}}.
#' @return Same structure as \code{\link{run_cmd}}: a named list with
#'   \code{status}, \code{stdout}, and \code{stderr}.
#' @export
run_conda_cmd <- function(conda_env,
                          cmd,
                          args = character(),
                          ...) {
  if (!nzchar(conda_env)) {
    rlang::abort(
      "conda_env must be a non-empty string.",
      class = "strainflow_invalid_arg"
    )
  }

  conda_base <- get_conda_base()

  # Build the full shell command string
  args_str <- if (length(args) > 0L) {
    paste(shQuote(args), collapse = " ")
  } else {
    ""
  }

  shell_cmd <- glue::glue(
    "source \"{conda_base}/etc/profile.d/conda.sh\" && ",
    "conda activate \"{conda_env}\" && ",
    "{cmd} {args_str}"
  )

  cli::cli_inform(c(
    "i" = "Running in conda env {.val {conda_env}}: {.code {cmd}}"
  ))

  run_cmd(
    cmd  = "bash",
    args = c("-c", shell_cmd),
    ...
  )
}

#' Get conda base path
#'
#' Determines the root directory of the conda installation by checking (in
#' order):
#' \enumerate{
#'   \item The \code{CONDA_EXE} environment variable (strips trailing
#'     \code{/bin/conda}).
#'   \item Common default installation paths on macOS and Linux.
#'   \item The \code{conda} executable found on \code{PATH}.
#' }
#'
#' @return Character string: absolute path to the conda base directory.
#' @noRd
get_conda_base <- function() {
  # 1. Check CONDA_EXE environment variable
  conda_exe <- Sys.getenv("CONDA_EXE", unset = "")
  if (nzchar(conda_exe) && fs::file_exists(conda_exe)) {
    # CONDA_EXE is typically /path/to/conda/bin/conda
    conda_base <- fs::path_dir(fs::path_dir(conda_exe))
    if (fs::file_exists(fs::path(conda_base, "etc", "profile.d", "conda.sh"))) {
      return(as.character(conda_base))
    }
  }

  # 2. Check common default locations
  home <- Sys.getenv("HOME", unset = "~")
  candidates <- c(
    fs::path(home, "miniconda3"),
    fs::path(home, "miniforge3"),
    fs::path(home, "anaconda3"),
    fs::path(home, "mambaforge"),
    "/opt/conda",
    "/opt/miniconda3",
    "/opt/homebrew/Caskroom/miniconda/base",
    "/usr/local/miniconda3"
  )

  for (candidate in candidates) {
    conda_sh <- fs::path(candidate, "etc", "profile.d", "conda.sh")
    if (fs::file_exists(conda_sh)) {
      return(as.character(candidate))
    }
  }

  # 3. Try to find conda on PATH and resolve from there
  conda_which <- Sys.which("conda")
  if (nzchar(conda_which)) {
    # Resolve symlinks
    resolved <- tryCatch(
      as.character(fs::path_real(conda_which)),
      error = function(e) conda_which
    )
    # conda binary is typically at <base>/bin/conda or <base>/condabin/conda
    candidate_base <- fs::path_dir(fs::path_dir(resolved))
    conda_sh <- fs::path(candidate_base, "etc", "profile.d", "conda.sh")
    if (fs::file_exists(conda_sh)) {
      return(as.character(candidate_base))
    }
    # Try one more level up (for condabin/conda)
    candidate_base2 <- fs::path_dir(candidate_base)
    conda_sh2 <- fs::path(candidate_base2, "etc", "profile.d", "conda.sh")
    if (fs::file_exists(conda_sh2)) {
      return(as.character(candidate_base2))
    }
  }

  rlang::abort(
    c(
      "Could not determine conda base directory.",
      "i" = "Set the CONDA_EXE environment variable to the path of your conda executable.",
      "i" = "Checked: CONDA_EXE env var, common default paths, and PATH."
    ),
    class = "strainflow_conda_base_not_found"
  )
}
