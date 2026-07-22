# Shared helpers for the portable NCI60 revision-code package.
#
# Copy this directory to the root of the NCI60 repository before running the
# companion scripts. The helpers deliberately do not refer to the manuscript
# revision workspace or to any other analysis repository.

revision_script_dir <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", arguments, value = TRUE)

  if (length(file_argument) > 0L) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_argument[[1]]),
      mustWork = TRUE
    )))
  }

  normalizePath(getwd(), mustWork = TRUE)
}

revision_repo_dir <- function(package_dir, env_name, required_paths) {
  override <- Sys.getenv(env_name, unset = "")
  candidates <- if (nzchar(override)) {
    override
  } else {
    c(dirname(package_dir), getwd())
  }

  for (candidate in candidates) {
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (dir.exists(candidate) && all(file.exists(file.path(candidate, required_paths)))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    "Could not locate the NCI60 repository. Copy revision_code_NCI60 to its ",
    "repository root or set ", env_name, ". Required paths: ",
    paste(required_paths, collapse = ", ")
  )
}

revision_output_dir <- function(package_dir) {
  override <- Sys.getenv("REVISION_OUTPUT_DIR", unset = "")
  output_dir <- if (nzchar(override)) override else file.path(package_dir, "outputs")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(output_dir, mustWork = TRUE)
}

revision_add_local_library <- function(repo_dir, env_name = "REVISION_R_LIB") {
  override <- Sys.getenv(env_name, unset = "")
  override_paths <- if (nzchar(override)) {
    strsplit(override, .Platform$path.sep, fixed = TRUE)[[1]]
  } else {
    character()
  }
  candidate_paths <- c(override_paths, file.path(repo_dir, ".r-lib"))
  existing_paths <- candidate_paths[dir.exists(candidate_paths)]

  if (length(existing_paths) > 0L) {
    .libPaths(unique(c(normalizePath(existing_paths), .libPaths())))
  }
}
