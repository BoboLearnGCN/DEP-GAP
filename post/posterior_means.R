#!/usr/bin/env Rscript

# Compute element-wise posterior means using MCMC files 5, 6, 7, 8, and 9.
#
# Expected input filenames in the directory supplied with -i include:
#   gamma_5.txt, ..., gamma_9.txt
#   mu_5.txt, ..., mu_9.txt
#   S_5.txt, ..., S_9.txt
#   sigma_square_5.txt, ..., sigma_square_9.txt
#
# Extensionless files and .csv files are also supported.
#
# Example:
#   Rscript compute_posterior_means.R -i MCMC_records -o posterior_results

parse_args <- function(args) {
  input_dir <- NULL
  output_dir <- NULL
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript compute_posterior_means.R -i INPUT_DIR -o OUTPUT_DIR\n\n",
        "Options:\n",
        "  -i, --input DIR    Directory containing the MCMC record files\n",
        "  -o, --output DIR   Directory for posterior-mean files\n",
        "  -h, --help         Show this help message\n",
        sep = ""
      )
      quit(status = 0L)
    } else if (arg %in% c("-i", "--input")) {
      if (i == length(args)) {
        stop("Missing directory after ", arg, call. = FALSE)
      }
      input_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (startsWith(arg, "--input=")) {
      input_dir <- sub("^--input=", "", arg)
      i <- i + 1L
    } else if (arg %in% c("-o", "--output")) {
      if (i == length(args)) {
        stop("Missing directory after ", arg, call. = FALSE)
      }
      output_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (startsWith(arg, "--output=")) {
      output_dir <- sub("^--output=", "", arg)
      i <- i + 1L
    } else {
      stop("Unknown argument: ", arg, "\nUse --help for usage.", call. = FALSE)
    }
  }

  if (is.null(input_dir) || !nzchar(input_dir)) {
    stop("The input directory is required. Example: -i MCMC_records", call. = FALSE)
  }
  if (!dir.exists(input_dir)) {
    stop("The input directory does not exist: ", input_dir, call. = FALSE)
  }
  if (is.null(output_dir) || !nzchar(output_dir)) {
    stop("The output directory is required. Example: -o posterior_results", call. = FALSE)
  }

  list(
    input_dir = normalizePath(input_dir, mustWork = TRUE),
    output_dir = normalizePath(output_dir, mustWork = FALSE)
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- args$input_dir
output_dir <- args$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

parameters <- c("gamma", "mu", "S", "sigma_square")
file_indices <- 5:9

find_input_file <- function(parameter, index) {
  stem <- file.path(input_dir, paste0(parameter, "_", index))
  candidates <- c(
    paste0(stem, ".txt"),
    paste0(stem, ".csv"),
    stem
  )
  existing <- candidates[file.exists(candidates)]

  if (!length(existing)) {
    stop(
      "Cannot find an input file for ", parameter, "_", index,
      ". Expected one of:\n  ", paste(candidates, collapse = "\n  "),
      call. = FALSE
    )
  }

  if (length(existing) > 1L) {
    stop(
      "Multiple input files found for ", parameter, "_", index, ":\n  ",
      paste(existing, collapse = "\n  "),
      "\nKeep only one version to avoid ambiguity.",
      call. = FALSE
    )
  }

  existing[[1L]]
}

read_numeric_file <- function(path) {
  extension <- tolower(tools::file_ext(path))

  if (extension == "csv") {
    value <- as.matrix(read.csv(path, header = FALSE, check.names = FALSE))
  } else {
    value <- as.matrix(read.table(
      path,
      header = FALSE,
      sep = "",
      fill = TRUE,
      check.names = FALSE
    ))
  }

  # Remove empty trailing columns that can be introduced by trailing delimiters.
  while (ncol(value) > 0L && all(is.na(value[, ncol(value)]))) {
    value <- value[, -ncol(value), drop = FALSE]
  }

  suppressWarnings(storage.mode(value) <- "double")
  if (!length(value) || anyNA(value)) {
    stop(
      "The file contains missing or non-numeric values: ", path,
      call. = FALSE
    )
  }

  value
}

posterior_mean <- function(parameter) {
  running_sum <- NULL
  expected_dim <- NULL

  for (index in file_indices) {
    path <- find_input_file(parameter, index)
    current <- read_numeric_file(path)

    if (is.null(running_sum)) {
      running_sum <- current
      expected_dim <- dim(current)
    } else {
      if (!identical(dim(current), expected_dim)) {
        stop(
          "Dimension mismatch for ", basename(path), ". Expected ",
          paste(expected_dim, collapse = " x "), "; found ",
          paste(dim(current), collapse = " x "), ".",
          call. = FALSE
        )
      }
      running_sum <- running_sum + current
    }
  }

  running_sum / length(file_indices)
}

message("Input directory:  ", input_dir)
message("Output directory: ", output_dir)
message("Using MCMC files with indices: ", paste(file_indices, collapse = ", "))

for (parameter in parameters) {
  message("Computing posterior mean for ", parameter, "...")
  result <- posterior_mean(parameter)
  destination <- file.path(output_dir, paste0(parameter, ".txt"))

  write.table(
    result,
    file = destination,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

message("Finished. All posterior means were saved in: ", output_dir)
