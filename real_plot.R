#!/usr/bin/env Rscript

# Combined from:
#   zhunbei_unormalized_unsupervised.Rmd
#   spatial_plots.Rmd
#
# Run from the directory containing the input data:
#   Rscript combined_analysis.R -o results

parse_args <- function(args) {
  output_dir <- NULL
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript combined_analysis.R -o OUTPUT_DIR\n",
        "\n",
        "Options:\n",
        "  -o, --output DIR   Directory for all generated files\n",
        "  -h, --help         Show this help message\n",
        sep = ""
      )
      quit(status = 0L)
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

  if (is.null(output_dir) || !nzchar(output_dir)) {
    stop("The output directory is required. Example: -o results", call. = FALSE)
  }

  list(output_dir = normalizePath(output_dir, mustWork = FALSE))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- normalizePath(getwd())
output_dir <- args$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "data.table", "dplyr", "ggplot2", "grid", "pheatmap", "reshape2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(pheatmap)
  library(reshape2)
})

input_file <- function(...) file.path(input_dir, ...)
output_file <- function(...) file.path(output_dir, ...)

assert_files_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Required input file(s) not found:\n  ",
      paste(missing, collapse = "\n  "),
      call. = FALSE
    )
  }
}

message("Input directory:  ", input_dir)
message("Output directory: ", output_dir)

# -----------------------------------------------------------------------------
# 1. Read posterior samples and calculate posterior summaries
# -----------------------------------------------------------------------------

K_max <- 10L
K_values <- 3:K_max
K_indices <- K_values - 2L
N_K <- length(K_values)

inversion_depth_file <- input_file("inversion_time_depth.txt")
neighbors_file <- input_file("neighbors_list.txt")
inversion_data_file <- input_file("zhunbei_cpp_2.0.txt")

parameters <- c("alpha", "mu", "gamma", "sigma_sq", "psi")
parameter_files <- lapply(parameters, function(parameter) {
  input_file(paste0("output_K", K_values), paste0(parameter, ".txt"))
})
names(parameter_files) <- parameters
post_label_files <- input_file(paste0("output_K", K_values), "post_labels.txt")

assert_files_exist(c(
  inversion_depth_file,
  neighbors_file,
  inversion_data_file,
  unlist(parameter_files, use.names = FALSE),
  post_label_files
))

inversion_depth <- read.table(inversion_depth_file, header = TRUE)
neighbors <- data.matrix(fread(neighbors_file, sep = " ", header = FALSE))
inversion_data <- array(
  unlist(fread(inversion_data_file, sep = " ", header = FALSE)),
  dim = c(3, 4600, 11)
)

record_list <- setNames(vector("list", length(parameters)), parameters)
post_list <- setNames(vector("list", length(parameters)), parameters)

for (parameter in parameters) {
  record_list[[parameter]] <- setNames(vector("list", N_K), paste0("K", K_values))
  post_list[[parameter]] <- setNames(vector("list", N_K), paste0("K", K_values))

  for (j in seq_along(K_values)) {
    samples <- data.matrix(fread(
      parameter_files[[parameter]][j], sep = " ", header = FALSE
    ))
    if (ncol(samples) < 5000L) {
      stop(
        basename(parameter_files[[parameter]][j]),
        " contains fewer than 5,000 recorded iterations.",
        call. = FALSE
      )
    }
    keep <- (ncol(samples) - 4999L):ncol(samples)
    record_list[[parameter]][[j]] <- samples
    post_list[[parameter]][[j]] <- rowMeans(samples[, keep, drop = FALSE])
  }
}

S_post <- setNames(vector("list", N_K), paste0("K", K_values))

for (j in seq_along(K_values)) {
  k <- K_values[j]
  post_list$mu[[j]] <- t(matrix(post_list$mu[[j]], nrow = k, ncol = 3))
  post_list$gamma[[j]] <- t(matrix(post_list$gamma[[j]], nrow = 11, ncol = 3))
  post_list$sigma_sq[[j]] <- t(matrix(
    post_list$sigma_sq[[j]], nrow = 11, ncol = 3
  ))

  psi <- matrix(0, nrow = k, ncol = k)
  index <- 1L
  for (k1 in seq_len(k - 1L)) {
    for (k2 in (k1 + 1L):k) {
      psi[k1, k2] <- psi[k2, k1] <- post_list$psi[[j]][index]
      index <- index + 1L
    }
  }
  post_list$psi[[j]] <- psi

  S_post[[j]] <- t(matrix(
    unlist(fread(post_label_files[j], sep = " ", header = FALSE)) + 1L,
    nrow = 4600,
    ncol = 11
  ))

  mu_summary <- post_list$mu[[j]]
  rownames(mu_summary) <- c("log(LLD)", "density", "velocity")
  colnames(mu_summary) <- paste0("lithology_", seq_len(k))
  write.csv(mu_summary, output_file(paste0("mu_K", k, ".csv")))

  gamma_summary <- post_list$gamma[[j]]
  rownames(gamma_summary) <- c("log(LLD)", "density", "velocity")
  colnames(gamma_summary) <- paste0("inversion_", 1:11)
  write.csv(gamma_summary, output_file(paste0("gamma_K", k, ".csv")))
}

# -----------------------------------------------------------------------------
# 2. Depth-effect plots and posterior heatmaps
# -----------------------------------------------------------------------------

depth_effects <- data.frame(
  data = unlist(post_list$gamma, use.names = FALSE),
  K = factor(rep(K_values, each = 33), levels = K_values),
  feature = rep(c("log_LLD", "density", "velocity"), times = 11L * N_K),
  depth = rep(rep(inversion_depth$depth, each = 3L), times = N_K)
)

feature_labels <- c(
  log_LLD = "log(LLD)", density = "density", velocity = "velocity"
)

for (feature_name in names(feature_labels)) {
  p <- ggplot(
    depth_effects[depth_effects$feature == feature_name, ],
    aes(x = depth, y = data, color = K, group = K)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    labs(
      y = feature_labels[[feature_name]],
      color = "# of lithology",
      title = paste("Depth effects on", feature_labels[[feature_name]])
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
    )

  ggsave(
    output_file(paste0("depth_effect_", feature_name, ".png")),
    p, width = 7, height = 5, dpi = 300, bg = "white"
  )
}

mu_k3 <- post_list$mu$K3
rownames(mu_k3) <- c("log(LLD)", "density", "velocity")
colnames(mu_k3) <- paste0("lithology_", 1:3)

pheatmap(
  mu_k3,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  cellwidth = 60,
  cellheight = 60,
  display_numbers = TRUE,
  fontsize_number = 15,
  number_format = "%.3f",
  angle_col = 45,
  fontsize_row = 12,
  fontsize_col = 12,
  filename = output_file("mu_K3_heatmap.png")
)

subtype_post_k6 <- sweep(post_list$mu$K6, 1, post_list$alpha$K6, FUN = "+")
norm_subtype_post_k6 <- t(scale(t(subtype_post_k6)))
rownames(norm_subtype_post_k6) <- c("log_LLD", "density", "velocity")
colnames(norm_subtype_post_k6) <- as.character(1:6)

pheatmap(
  norm_subtype_post_k6,
  cluster_rows = FALSE,
  filename = output_file("normalized_subtype_K6_heatmap.png")
)

physical_inversion <- data.frame(
  data = c(t(norm_subtype_post_k6)),
  type = rep(rownames(norm_subtype_post_k6), each = 6),
  lithology = factor(rep(1:6, 3), levels = c(1, 4, 2, 6, 3, 5))
)

p <- ggplot(
  physical_inversion,
  aes(x = lithology, y = data, group = type, color = type)
) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "predicted lithology", y = "normalized values", color = "feature") +
  theme_bw()

ggsave(
  output_file("normalized_physical_inversion_K6.png"),
  p, width = 7, height = 5, dpi = 300, bg = "white"
)

# -----------------------------------------------------------------------------
# 3. BIC and BLIC
# -----------------------------------------------------------------------------

BIC_mixture <- function(data_obs, .K, .B, .n, .alpha, .mu, .gamma,
                        .sigma_sq, .S) {
  value <- 0
  log_terms <- matrix(NA_real_, nrow = .n, ncol = .K)
  proportions_by_layer <- t(apply(.S, 1, function(labels) {
    tab <- table(factor(labels, levels = seq_len(.K)))
    as.numeric(tab) / sum(tab)
  }))

  for (b in seq_len(.B)) {
    for (k in seq_len(.K)) {
      mean_bk <- .alpha + .mu[, k] + .gamma[, b]
      sd_bk <- sqrt(.sigma_sq[, b])
      log_terms[, k] <- log(proportions_by_layer[b, k]) +
        colSums(dnorm(data_obs[, , b], mean = mean_bk, sd = sd_bk, log = TRUE))
    }
    row_max <- apply(log_terms, 1, max)
    value <- value + sum(row_max) + sum(log(rowSums(exp(log_terms - row_max))))
  }
  value
}

BLIC <- function(data_obs, .K, .B, .n, .alpha, .mu, .gamma, .sigma_sq,
                 .psi, .S, .neighbors) {
  value <- 0
  log_terms <- matrix(NA_real_, nrow = .n, ncol = .K)
  potts_terms <- matrix(NA_real_, nrow = .n, ncol = .K)
  for (b in seq_len(.B)) {
    labels <- as.numeric(.S[b, ])
    for (k in seq_len(.K)) {
      mean_bk <- .alpha + .mu[, k] + .gamma[, b]
      sd_bk <- sqrt(.sigma_sq[, b])
      potts_terms[, k] <- apply(.neighbors, 2, function(idx) {
        valid_idx <- idx[idx > 0]
        sum(.psi[labels[valid_idx], k])
      })
      log_terms[, k] <- colSums(
        dnorm(data_obs[, , b], mean = mean_bk, sd = sd_bk, log = TRUE)
      )
    }

    potts_max <- apply(potts_terms, 1, max)
    potts_prob <- exp(potts_terms - potts_max)
    potts_prob <- potts_prob / rowSums(potts_prob)
    log_terms <- log_terms + log(potts_prob)
    row_max <- apply(log_terms, 1, max)
    value <- value + sum(row_max) + sum(log(rowSums(exp(log_terms - row_max))))
  }
  value
}

G <- 3L
B <- 11L
n <- 4600L
BIC <- numeric(N_K)
BLIC_values <- numeric(N_K)

for (j in seq_along(K_values)) {
  k <- K_values[j]
  parameter_count <- k * G + (2 * B - 1) * G + k * (k - 1) / 2

  BIC[j] <- -2 * BIC_mixture(
    inversion_data, k, B, n,
    post_list$alpha[[j]], post_list$mu[[j]], post_list$gamma[[j]],
    post_list$sigma_sq[[j]], S_post[[j]]
  ) + parameter_count * log(G * B * n)

  BLIC_values[j] <- -2 * BLIC(
    inversion_data, k, B, n,
    post_list$alpha[[j]], post_list$mu[[j]], post_list$gamma[[j]],
    post_list$sigma_sq[[j]], post_list$psi[[j]], S_post[[j]], neighbors
  ) + parameter_count * log(G * B * n)

  message("K = ", k, ": BIC = ", BIC[j], "; BLIC = ", BLIC_values[j])
}

criteria <- data.frame(K = K_values, BIC = BIC, BLIC = BLIC_values)
write.csv(criteria, output_file("information_criteria.csv"), row.names = FALSE)

criteria_long <- reshape2::melt(
  criteria, id.vars = "K", variable.name = "criterion", value.name = "value"
)
p <- ggplot(criteria_long, aes(x = K, y = value, color = criterion)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  labs(x = "K", y = "Information criterion", color = NULL) +
  theme_bw() +
  theme(axis.title = element_text(size = 15))

ggsave(
  output_file("BIC_BLIC.png"),
  p, width = 7, height = 5, dpi = 300, bg = "white"
)

# -----------------------------------------------------------------------------
# 4. Spatial label maps for K = 6
# -----------------------------------------------------------------------------

coordinates_file <- input_file("inversion_coordinates.txt")
spatial_labels_file <- input_file("post_labels_K6.txt")
assert_files_exist(c(coordinates_file, spatial_labels_file))

coord <- read.table(coordinates_file, header = FALSE)
post_labels_k6 <- data.matrix(fread(spatial_labels_file, sep = " ", header = FALSE))

color_box <- c("#50BF50", "#E8585D", "#FFC0CB", "#6495ED", "#FFFF00", "#FFA500")
x_range <- diff(range(coord[, 1], na.rm = TRUE))
y_range <- diff(range(coord[, 2], na.rm = TRUE))

for (d in seq_len(nrow(post_labels_k6))) {
  plot_data <- data.frame(
    x = coord[, 1],
    y = coord[, 2],
    label = factor(unlist(post_labels_k6[d, ]))
  )

  p <- ggplot(plot_data, aes(x = x, y = y, color = label)) +
    geom_point(show.legend = FALSE, size = 0.8, shape = 15) +
    scale_color_manual(values = color_box) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed(expand = FALSE) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "none"
    )

  ggsave(
    output_file(paste0("spatial_map_", d, ".png")),
    p,
    width = 1.2 * x_range / y_range,
    height = 1.2,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

# -----------------------------------------------------------------------------
# 5. Well-log plot
# -----------------------------------------------------------------------------

well_file <- normalizePath(
  file.path(input_dir, "..", "well_data", "zhunbei_processed.csv"),
  mustWork = FALSE
)
assert_files_exist(well_file)

zhunbei_well <- read.csv(well_file)
ZB3 <- zhunbei_well[zhunbei_well$Well == "ZB3", ]
ZB3$lithology_level2_correct <- ifelse(
  ZB3$lithology_level2 == "TYY", "CM", "other"
)

p <- ggplot(
  ZB3,
  aes(
    x = log_LLD,
    y = DEPTH / 1000,
    color = lithology_level2_correct,
    size = lithology_level2_correct,
    alpha = lithology_level2_correct
  )
) +
  geom_point() +
  scale_y_reverse() +
  scale_x_continuous(position = "top") +
  scale_color_manual(values = c(CM = "#FFFF00", other = "#CFCFCF")) +
  scale_size_manual(values = c(CM = 1, other = 0.5)) +
  scale_alpha_manual(values = c(CM = 1, other = 0.5)) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 32),
    axis.text = element_text(size = 26)
  ) +
  labs(x = "log(LLD)", y = "Depth (km)")

ggsave(
  output_file("depth_logLLD_points.png"),
  p, width = 5, height = 8, dpi = 300, bg = "white"
)

# -----------------------------------------------------------------------------
# 6. Revised lithology-effect heatmap from spatial_plots.Rmd
# -----------------------------------------------------------------------------

mat <- matrix(
  c(
    1.45, 1.55, 1.85, 1.25, 1.60, 1.35,
    2.00, 1.25, 1.95, 1.95, 1.25, 0.95,
    0.65, 1.10, 0.85, 0.35, 0.00, 0.20
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(paste0("R", 1:3), c("red", "blue", "green", "orange", "pink", "yellow"))
)

annotation <- data.frame(
  group = c("red", "blue", "green", "orange", "pink", "yellow"),
  col = c("#E8585D", "#6495ED", "#50BF50", "#FFA500", "#FFC0CB", "#FFFF00")
)
new_order <- c("green", "red", "pink", "blue", "orange", "yellow")
mat <- mat[, new_order]
annotation <- annotation[match(new_order, annotation$group), ]
annotation$x <- seq_len(nrow(annotation))
mat <- sweep(mat, 1, mat[, 1], FUN = "-")
mat <- pmax(-1, pmin(1, mat))

heat_df <- as.data.frame(as.table(mat))
colnames(heat_df) <- c("Row", "Col", "Value")
heat_df$x <- as.numeric(heat_df$Col)
heat_df$y <- 4 - as.numeric(heat_df$Row)

p <- ggplot() +
  geom_tile(
    data = heat_df,
    aes(x = x, y = y, fill = Value),
    width = 1, height = 1, color = NA
  ) +
  geom_rect(
    data = annotation,
    aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 3.58, ymax = 3.72),
    fill = annotation$col,
    color = NA,
    inherit.aes = FALSE
  ) +
  annotate(
    "rect", xmin = 0.5, xmax = 6.5, ymin = 3.50, ymax = 3.58,
    fill = "white", color = NA
  ) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-1, 1), breaks = c(-1, 0, 1),
    guide = guide_colorbar(
      title = NULL, direction = "horizontal",
      barwidth = unit(1.5, "cm"), barheight = unit(0.25, "cm"),
      ticks = TRUE, label.position = "bottom"
    )
  ) +
  coord_fixed(clip = "off") +
  scale_x_continuous(limits = c(0.5, 6.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 4.45), expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = c(0.80, 0.91),
    legend.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(
  output_file("mu_heatmap_revised.png"),
  p, width = 5.5, height = 3.5, dpi = 300, bg = "white"
)

# -----------------------------------------------------------------------------
# 7. Compact depth-effect heatmap from spatial_plots.Rmd
# -----------------------------------------------------------------------------

depth_mat <- matrix(
  c(
    -0.07, 0.03, 0.08, 0.18, 0.24, 0.24, 0.21, 0.24, 0.24, 0.29, 0.29,
    -0.07, -0.16, -0.22, -0.16, -0.02, 0.03, 0.08, -0.07, -0.16, -0.26, -0.30,
    -0.07, -0.19, -0.07, 0.08, 0.30, 0.30, -0.19, -0.24, -0.30, -0.07, 0.18
  ),
  nrow = 3,
  byrow = TRUE
)

depth_df <- as.data.frame(as.table(depth_mat))
colnames(depth_df) <- c("Row", "Col", "Value")
depth_df$x <- as.numeric(depth_df$Col)
depth_df$y <- nrow(depth_mat) - as.numeric(depth_df$Row) + 1

depth_annotation <- data.frame(
  x = seq_len(ncol(depth_mat)),
  col = colorRampPalette(c("#00C8D7", "#3B95D9", "#7B65E8", "#C332F4", "#FF00C8"))(
    ncol(depth_mat)
  )
)

p <- ggplot() +
  geom_tile(
    data = depth_df,
    aes(x = x, y = y, fill = Value),
    width = 1, height = 1, color = NA
  ) +
  geom_rect(
    data = depth_annotation,
    aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 3.56, ymax = 3.70),
    fill = depth_annotation$col,
    color = NA,
    inherit.aes = FALSE
  ) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_fixed(clip = "off") +
  scale_x_continuous(limits = c(0.5, ncol(depth_mat) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 3.72), expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave(
  output_file("depth_effects_heatmap.png"),
  p, width = 1.5, height = 0.5, dpi = 300, bg = "white"
)

message("Analysis complete. Results were written to: ", output_dir)
