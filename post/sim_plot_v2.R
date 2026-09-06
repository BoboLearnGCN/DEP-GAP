#!/usr/bin/env Rscript

# Posterior analysis for K = 6, B = 9, n = 100000 per layer
# Run from the directory containing alpha.txt, mu.txt, gamma.txt, etc.:
#   Rscript unsupervised_K6_revised.R -o results

parse_args <- function(args) {
  output_dir <- NULL
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript unsupervised_K6_revised.R -o OUTPUT_DIR\n\n",
        "Options:\n",
        "  -o, --output DIR   Directory for every generated file\n",
        "  -h, --help         Show this help message\n",
        sep = ""
      )
      quit(status = 0L)
    } else if (arg %in% c("-o", "--output")) {
      if (i == length(args) || startsWith(args[[i + 1L]], "-")) {
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
if (!dir.exists(output_dir)) {
  stop("Could not create output directory: ", output_dir, call. = FALSE)
}

input_file <- function(...) file.path(input_dir, ...)
output_file <- function(...) file.path(output_dir, ...)

required_packages <- c(
  "agricolae", "data.table", "ggplot2", "gridExtra", "heatmap3",
  "pheatmap", "reshape2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Install these R packages before running the script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(agricolae)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(heatmap3)
  library(pheatmap)
  library(reshape2)
})

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

read_numeric_matrix <- function(path) {
  # Read wide, whitespace-delimited matrices without creating one column object
  # for each observation. Rows in each file remain rows in the returned matrix.
  fields <- count.fields(
    path, sep = "", quote = "", comment.char = "", blank.lines.skip = TRUE
  )
  if (!length(fields) || anyNA(fields) || any(fields < 1L) ||
      any(fields != fields[1L])) {
    stop(basename(path), " must be a rectangular numeric matrix.", call. = FALSE)
  }
  values <- scan(path, what = double(), quiet = TRUE, comment.char = "")
  if (length(values) != sum(fields) || any(!is.finite(values))) {
    stop(basename(path), " contains missing, non-finite, or malformed values.",
         call. = FALSE)
  }
  matrix(values, nrow = length(fields), ncol = fields[1L], byrow = TRUE)
}

assert_dimensions <- function(x, expected, name) {
  if (!identical(dim(x), as.integer(expected))) {
    stop(name, " has dimensions ", paste(dim(x), collapse = " x "),
         "; expected ", paste(expected, collapse = " x "), ".",
         call. = FALSE)
  }
}

read_label_matrix <- function(path, expected_rows, expected_columns, K) {
  labels <- read_numeric_matrix(path)
  assert_dimensions(labels, c(expected_rows, expected_columns), basename(path))
  if (any(labels != floor(labels))) {
    stop(basename(path), " must contain integer labels.", call. = FALSE)
  }
  if (any(labels == 0)) {
    labels <- labels + 1L
  }
  if (any(!labels %in% seq_len(K))) {
    stop(basename(path), " labels must be between 1 and ", K,
         " after conversion.", call. = FALSE)
  }
  storage.mode(labels) <- "integer"
  labels
}

message("Input directory:  ", input_dir)
message("Output directory: ", output_dir)

# -----------------------------------------------------------------------------
# 1. Load posterior records and input data
# -----------------------------------------------------------------------------

alpha_file <- input_file("alpha.txt")
mu_file <- input_file("mu.txt")
gamma_file <- input_file("gamma.txt")
sigma_sq_file <- input_file("sigma_sq.txt")
psi_file <- input_file("psi.txt")
observed_file <- input_file("Y_obs.txt")
labels_file <- input_file("Z.txt")
coordinates_file <- input_file("coords.csv")

assert_files_exist(c(
  alpha_file, mu_file, gamma_file, sigma_sq_file, psi_file,
  observed_file, labels_file, coordinates_file
))

K <- 6L
B <- 9L
n <- 100000L
n_iter <- 10000L
n_burnin <- 5000L
after_burnin <- seq.int(n_burnin + 1L, n_iter)
n_record <- length(after_burnin)

alpha_record <- read_numeric_matrix(alpha_file)
mu_record <- read_numeric_matrix(mu_file)
gamma_record <- read_numeric_matrix(gamma_file)
sigma_sq_record <- read_numeric_matrix(sigma_sq_file)
psi_record <- read_numeric_matrix(psi_file)
Y_obs <- read_numeric_matrix(observed_file)
S_post <- read_label_matrix(labels_file, B, n, K)
G <- nrow(alpha_record)

inversion_coordinates <- as.data.frame(fread(coordinates_file, header = "auto"))
if (ncol(inversion_coordinates) != 2L || nrow(inversion_coordinates) != n) {
  stop("coords.csv must contain ", n, " rows and exactly two columns.",
       call. = FALSE)
}
names(inversion_coordinates) <- c("x", "y")
if (!all(vapply(inversion_coordinates, is.numeric, logical(1))) ||
    any(!is.finite(as.matrix(inversion_coordinates)))) {
  stop("coords.csv must contain finite numeric coordinates.", call. = FALSE)
}

property_names <- paste0("P_", seq_len(G))
rocktype_names <- paste0("rocktype", seq_len(K))
layer_names <- paste0("L", seq_len(B))

records <- list(
  alpha = alpha_record,
  mu = mu_record,
  gamma = gamma_record,
  sigma_sq = sigma_sq_record,
  psi = psi_record
)
too_short <- names(records)[vapply(records, ncol, integer(1)) < n_iter]
if (length(too_short)) {
  stop(
    "The following MCMC files contain fewer than 10,000 columns: ",
    paste(too_short, collapse = ", "),
    call. = FALSE
  )
}

assert_dimensions(Y_obs, c(G, n * B), "Y_obs.txt")
assert_dimensions(S_post, c(B, n), "Z.txt")
expected_rows <- c(alpha = G, mu = G * K, gamma = G * B,
                   sigma_sq = G * B, psi = choose(K, 2L))
for (record_name in names(records)) {
  if (nrow(records[[record_name]]) != expected_rows[[record_name]]) {
    stop(record_name, ".txt has ", nrow(records[[record_name]]),
         " rows; expected ", expected_rows[[record_name]], ".", call. = FALSE)
  }
}

rownames(Y_obs) <- property_names
colnames(Y_obs) <- paste0(rep(seq_len(B), each = n), "_", rep(seq_len(n), B))
rownames(S_post) <- layer_names

# Y_obs is G x (n * B); each element of Y_obs_list is G x n.
Y_obs_list <- lapply(seq_len(B), function(b) {
  idx <- (b - 1L) * n + seq_len(n)
  Y_obs[, idx, drop = FALSE]
})
names(Y_obs_list) <- layer_names

Y_melt <- reshape2::melt(Y_obs)
colnames(Y_melt) <- c("property", "sample", "value")
Y_melt$layer <- factor(rep(layer_names, each = n * G), levels = layer_names)

# -----------------------------------------------------------------------------
# 2. Posterior summaries and summary CSV files
# -----------------------------------------------------------------------------

alpha_post <- rowMeans(alpha_record[, after_burnin, drop = FALSE])
# Records are grouped by property, then rock type or layer.
mu_post <- t(matrix(
  rowMeans(mu_record[, after_burnin, drop = FALSE]), nrow = K, ncol = G
))
gamma_post <- t(matrix(
  rowMeans(gamma_record[, after_burnin, drop = FALSE]), nrow = B, ncol = G
))
sigma_sq_post <- t(matrix(
  rowMeans(sigma_sq_record[, after_burnin, drop = FALSE]), nrow = B, ncol = G
))
rocktype_post <- sweep(mu_post, 1, alpha_post, FUN = "+")

names(alpha_post) <- property_names
rownames(mu_post) <- property_names
colnames(mu_post) <- rocktype_names
rownames(rocktype_post) <- property_names
colnames(rocktype_post) <- rocktype_names
rownames(gamma_post) <- property_names
colnames(gamma_post) <- layer_names
rownames(sigma_sq_post) <- property_names
colnames(sigma_sq_post) <- layer_names
assert_dimensions(mu_post, c(G, K), "mu_post")
assert_dimensions(gamma_post, c(G, B), "gamma_post")
if (any(!is.finite(sigma_sq_post)) || any(sigma_sq_post <= 0)) {
  stop("Posterior variances must be finite and positive.", call. = FALSE)
}

write.csv(alpha_post, output_file("alpha_K6.csv"))
write.csv(mu_post, output_file("mu_K6.csv"))
write.csv(rocktype_post, output_file("rocktype_K6.csv"))
write.csv(gamma_post, output_file("gamma_K6.csv"))
write.csv(sigma_sq_post, output_file("sigma_sq_K6.csv"))

# -----------------------------------------------------------------------------
# 3. Parameter heatmaps
# -----------------------------------------------------------------------------

rocktype_est <- mu_post
rownames(rocktype_est) <- property_names
colnames(rocktype_est) <- rocktype_names
rocktype_melt <- reshape2::melt(rocktype_est)
colnames(rocktype_melt) <- c("property", "rocktype", "value")
rocktype_melt$rocktype <- factor(rocktype_melt$rocktype, levels = rocktype_names)

p <- ggplot(rocktype_melt, aes(x = rocktype, y = property, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3))) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-1.5, 1.5), oob = scales::squish
  ) +
  scale_x_discrete(position = "top") +
  labs(title = "Estimated rock type effects, K=6", x = "rock type", y = "property") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
ggsave(output_file("estimated_rocktype_effects.png"), p, width = 7, height = 4, dpi = 300)

gamma_plot_matrix <- gamma_post
rownames(gamma_plot_matrix) <- property_names
colnames(gamma_plot_matrix) <- layer_names
gamma_melt <- reshape2::melt(gamma_plot_matrix)
colnames(gamma_melt) <- c("property", "layer", "value")
gamma_melt$layer <- factor(gamma_melt$layer, levels = colnames(gamma_plot_matrix))

p <- ggplot(gamma_melt, aes(x = layer, y = property, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3))) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-0.35, 0.35), oob = scales::squish
  ) +
  scale_x_discrete(position = "top") +
  labs(title = "Estimated layer effects", x = "layer", y = "property") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
ggsave(output_file("estimated_layer_effects.png"), p, width = 10, height = 4, dpi = 300)

sigma_melt <- reshape2::melt(sigma_sq_post)
colnames(sigma_melt) <- c("property", "layer", "value")
sigma_melt$layer <- factor(sigma_melt$layer, levels = layer_names)
p <- ggplot(sigma_melt, aes(x = layer, y = property, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  scale_x_discrete(position = "top") +
  labs(title = "Estimated variance effects", x = "layer", y = "property") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
ggsave(output_file("estimated_layer_variance_effects.png"), p, width = 9, height = 4, dpi = 300)

# Publication-style parameter tables from the notebook
rocktype_app_df <- reshape2::melt(rocktype_post)
colnames(rocktype_app_df) <- c("property", "rocktype", "value")
rocktype_app_df$rocktype <- factor(rocktype_app_df$rocktype, levels = rocktype_names)
p <- ggplot(rocktype_app_df, aes(x = rocktype, y = property, fill = value)) +
  geom_tile() +
  geom_text(aes(label = format(round(value, 2), nsmall = 2)), size = 7, color = "green3") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "", y = "") +
  theme(
    axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 22),
    legend.text = element_text(size = 20), legend.title = element_text(size = 22)
  )
ggsave(output_file("rocktype_app.png"), p, dpi = 100, width = 6, height = 3, units = "in")

gamma_app_df <- gamma_melt
p <- ggplot(gamma_app_df, aes(x = layer, y = property, fill = value)) +
  geom_tile() +
  geom_text(aes(label = format(round(value, 2), nsmall = 2)), size = 6, color = "green3") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(
    low = "skyblue2", mid = "white", high = "#E3095C",
    midpoint = 0, limits = c(-0.35, 0.35), oob = scales::squish
  ) +
  labs(x = "", y = "") +
  theme(
    axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 20),
    legend.text = element_text(size = 20), legend.title = element_text(size = 22)
  )
ggsave(output_file("gamma_app.png"), p, dpi = 100, width = 9, height = 3, units = "in")

# heatmap3 versions
if (G >= 2L) {
  rocktype_heatmap <- rocktype_post[rev(seq_len(G)), , drop = FALSE]
  colnames(rocktype_heatmap) <- paste0("type ", seq_len(K))
  png(output_file("est_rocktype.png"), width = 600, height = 500)
  heatmap3(
    rocktype_heatmap, Colv = NA, Rowv = NA, scale = "none",
    ColSideColors = c("#E8585D", "cornflowerblue", "#50BF50", "orange", "pink", "yellow"),
    col = colorRampPalette(c("white", "black"))(20),
    ColSideLabs = FALSE, labRow = NA, labCol = NA
  )
  dev.off()

  layer_heatmap <- gamma_plot_matrix[rev(seq_len(G)), , drop = FALSE]
  colnames(layer_heatmap) <- paste0("L", seq_len(B))
  gradient_colors <- colorRampPalette(c("#00dbde", "#fc00ff"))(B)
  png(output_file("est_layer.png"), width = 600, height = 500)
  heatmap3(
    layer_heatmap, Colv = NA, Rowv = NA, scale = "none",
    ColSideColors = gradient_colors,
    col = colorRampPalette(c("blue", "white", "red"))(20),
    ColSideLabs = FALSE, labRow = NA, labCol = NA
  )
  dev.off()
}

# -----------------------------------------------------------------------------
# 4. Spatial clustering maps
# -----------------------------------------------------------------------------

color_box <- c("#E8585D", "cornflowerblue", "#50BF50", "orange", "pink", "yellow")

for (b in seq_len(B)) {
  spatial_df <- data.frame(
    x = inversion_coordinates$x,
    y = inversion_coordinates$y,
    label = factor(as.numeric(S_post[b, ]), levels = seq_len(K))
  )

  p <- ggplot(spatial_df, aes(x = x, y = y, color = label)) +
    geom_point(size = 0.8, shape = 15, show.legend = FALSE) +
    scale_color_manual(values = color_box, drop = FALSE) +
    coord_fixed() +
    labs(x = "latitude", y = "longitude", color = "label") +
    theme_bw() +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.title = element_text(size = 20), panel.border = element_blank(),
      panel.grid = element_blank()
    )

  ggsave(
    output_file(paste0("Clustering_K6_", b, ".png")),
    p, units = "in", width = 5, height = 2.5, dpi = 150
  )
}

# -----------------------------------------------------------------------------
# 5. Observed-and-corrected-data plots by layer and predicted label
# -----------------------------------------------------------------------------

S_vec <- as.vector(t(S_post))
Y_melt$label <- factor(rep(S_vec, each = G), levels = seq_len(K))

p <- ggplot(
  Y_melt[Y_melt$property == "P_1", ],
  aes(x = label, y = value, fill = factor(layer))
) +
  geom_boxplot() +
  labs(x = "predicted label", y = "P1", fill = "layer") +
  theme_bw()
ggsave(output_file("P1_by_label_and_layer.png"), p, width = 10, height = 6, dpi = 300)

gradient_colors <- colorRampPalette(c("#00dbde", "#fc00ff"))(B)
label1_resistivity <- Y_melt[
  Y_melt$property == "P_1" & Y_melt$label == "1", ]
label1_median <- if (nrow(label1_resistivity)) {
  aggregate(value ~ layer, data = label1_resistivity, FUN = median)
} else {
  data.frame(layer = character(), value = numeric())
}
layer_levels <- rev(layer_names)
label1_resistivity$layer <- factor(label1_resistivity$layer, levels = layer_levels)
label1_median$layer <- factor(label1_median$layer, levels = layer_levels)

p <- ggplot() +
  geom_boxplot(
    data = label1_resistivity,
    aes(x = layer, y = value, fill = layer), show.legend = FALSE
  ) +
  geom_line(
    data = label1_median,
    aes(x = layer, y = value, group = 1),
    color = "#FFD700", linewidth = 1.2, linetype = "longdash"
  ) +
  scale_fill_manual(values = gradient_colors) +
  scale_y_continuous(position = "right") +
  coord_flip() +
  labs(x = "layer", y = "P_1") +
  theme_bw() +
  theme(
    panel.grid = element_blank(), panel.border = element_blank(),
    axis.line = element_line(color = "black"), axis.title = element_text(size = 20),
    axis.text = element_text(size = 14)
  )
ggsave(output_file("P1_label1_by_layer.png"), p, width = 7, height = 6, dpi = 300)

# Scatter plots of observed properties
scatter_labels <- factor(as.vector(t(S_post)), levels = seq_len(K))
if (G >= 2L) {
  p <- ggplot(data.frame(R = Y_obs[1, ], D = Y_obs[2, ], label = scatter_labels),
              aes(x = R, y = D, color = label)) +
    geom_point(size = 0.4) +
    labs(color = "label", y = "P2", x = "P1", title = "I1")
  ggsave(output_file("P1_vs_P2.png"), p, width = 7, height = 5, dpi = 300)
}

if (G >= 3L) {
  p <- ggplot(data.frame(R = Y_obs[1, ], V = Y_obs[3, ], label = scatter_labels),
              aes(x = R, y = V, color = label)) +
    geom_point(size = 0.4) +
    labs(color = "label", y = "P3", x = "P1", title = "I1")
  ggsave(output_file("P1_vs_P3.png"), p, width = 7, height = 5, dpi = 300)
}

spatial_obs <- vector("list", B)
for (b in seq_len(B)) {
  spatial_obs[[b]] <- data.frame(
    x = inversion_coordinates$x,
    y = inversion_coordinates$y,
    value = Y_obs_list[[b]][1, ]
  )
}

spatial_heat_obs <- vector("list", B)

for (b in seq_len(B)) {
  spatial_heat_obs[[b]] <- ggplot(spatial_obs[[b]], aes(x=x, y=y, fill=value)) +
    geom_tile(show.legend = FALSE) +
    scale_fill_gradient(low = "white", high = "black", limits=c(4.4, 6.7), oob = scales::squish) +
    coord_fixed() +
    theme_bw() +
    theme(axis.text = element_blank(), axis.title = element_text(size = 36), legend.title = element_text(size = 22, vjust = 0.1), legend.text = element_text(size = 20), panel.grid = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  
  ggsave(filename = output_file(paste0("spatial_heat_obs_", b, ".png")), spatial_heat_obs[[b]], dpi = 100, height = 7, width = 7)
}

Y_correct_list <- vector("list", B)
names(Y_correct_list) <- layer_names
spatial_correct <- vector("list", B)
for (b in seq_len(B)) {
  # Center on alpha + mu, remove the layer mean effect, and match layer 1 variance.
  baseline <- rocktype_post[, S_post[b, ], drop = FALSE]
  residual <- sweep(Y_obs_list[[b]] - baseline, 1L, gamma_post[, b], "-")
  multiplier <- sqrt(sigma_sq_post[, 1L] / sigma_sq_post[, b])
  Y_correct_list[[b]] <- baseline + sweep(residual, 1L, multiplier, "*")
  dimnames(Y_correct_list[[b]]) <- dimnames(Y_obs_list[[b]])
  spatial_correct[[b]] <- data.frame(
    x = inversion_coordinates$x,
    y = inversion_coordinates$y,
    value = Y_correct_list[[b]][1, ]
  )
}

spatial_heat_correct <- vector("list", B)

for (b in seq_len(B)) {
  spatial_heat_correct[[b]] <-
    ggplot(spatial_correct[[b]],aes(x=x, y=y, fill=value)) +
    geom_tile(show.legend = FALSE) +
    scale_fill_gradient(low = "white", high = "black", limits=c(4.4, 6.7), oob = scales::squish) +
    coord_fixed()+
    theme_bw() +
    theme(axis.text = element_blank(), axis.title = element_text(size = 36, hjust = 0.5, vjust = 0.5), legend.title = element_text(size = 22, vjust = 0.1), legend.text = element_text(size = 20), panel.grid = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  
  ggsave(filename = output_file(paste0("spatial_heat_correct_", b, ".png")), spatial_heat_correct[[b]], dpi = 100, height = 7, width = 7)
}


# -----------------------------------------------------------------------------
# 6. Psi MCMC diagnostics, tests, and credible intervals
# -----------------------------------------------------------------------------

n_psi <- choose(K, 2)
if (nrow(psi_record) != n_psi) {
  stop("psi.txt must contain 15 rows for K=6; found ", nrow(psi_record), ".", call. = FALSE)
}

psi_index <- c(
  paste0("psi_1", 2:6), paste0("psi_2", 3:6), paste0("psi_3", 4:6),
  paste0("psi_4", 5:6), "psi_56"
)
psi_sample <- psi_record[, after_burnin, drop = FALSE]
psi_post <- rowMeans(psi_sample)
# Fraction of saved post-burn-in transitions that change the recorded value.
# This equals acceptance only when every proposal is recorded without rounding.
acceptance_rate <- rowMeans(
  psi_record[, after_burnin, drop = FALSE] !=
    psi_record[, after_burnin - 1L, drop = FALSE]
)
write.csv(
  data.frame(parameter = psi_index, posterior_mean = psi_post, acceptance_rate = acceptance_rate),
  output_file("psi_summary.csv"), row.names = FALSE
)

psi_plots <- lapply(seq_len(n_psi), function(i) {
  trace_df <- data.frame(iteration = after_burnin, psi = psi_sample[i, ])
  ggplot(trace_df, aes(x = iteration, y = psi)) +
    geom_line(linewidth = 0.3) +
    geom_hline(yintercept = mean(trace_df$psi), color = "green") +
    labs(title = paste("MCMC plot of", psi_index[i]), x = "iteration", y = psi_index[i]) +
    theme_bw() +
    theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))
})
trace_grob <- gridExtra::arrangeGrob(grobs = psi_plots, ncol = 1)
ggsave(output_file("psi_MCMC_traces.png"), trace_grob, width = 8, height = 30, dpi = 200, limitsize = FALSE)

# Legacy exploratory summaries retained from the original analysis.
# Posterior draws are not independent observed replicates: these t/normal/LSD
# p-values are not calibrated significance tests for differences in psi.
# Use the paired posterior intervals below for posterior uncertainty.
psi_combn <- combn(n_psi, 2)
n_comparisons <- ncol(psi_combn)
mean_diff <- numeric(n_comparisons)
sd_pool <- numeric(n_comparisons)
psi_test_stat <- numeric(n_comparisons)
psi_test_stat_2 <- numeric(n_comparisons)

for (i in seq_len(n_comparisons)) {
  idx1 <- psi_combn[1, i]
  idx2 <- psi_combn[2, i]
  mean_diff[i] <- psi_post[idx1] - psi_post[idx2]
  sd_pool[i] <- sqrt(
    (var(psi_sample[idx1, ]) * (n_record - 1) +
       var(psi_sample[idx2, ]) * (n_record - 1)) /
      (2 * n_record - 2)
  )
  denominator <- sd_pool[i] * sqrt(2 / n_record)
  psi_test_stat[i] <- if (denominator > 0) mean_diff[i] / denominator else
    if (mean_diff[i] == 0) 0 else sign(mean_diff[i]) * Inf
  delta_sd <- sd(psi_sample[idx1, ] - psi_sample[idx2, ])
  psi_test_stat_2[i] <- if (delta_sd > 0) mean_diff[i] / delta_sd else
    if (mean_diff[i] == 0) 0 else sign(mean_diff[i]) * Inf
}

psi_p_value <- 2 * pt(-abs(psi_test_stat), df = 2 * n_record - 2)
psi_p_value_2 <- p.adjust(2 * pnorm(-abs(psi_test_stat_2)), method = "BH")
comparison_labels <- apply(psi_combn, 2, function(idx) {
  paste(psi_index[idx[1]], "vs", psi_index[idx[2]])
})
psi_tests <- data.frame(
  comparison = comparison_labels,
  mean_difference = mean_diff,
  pooled_sd = sd_pool,
  t_statistic = psi_test_stat,
  t_p_value = psi_p_value,
  standardized_statistic = psi_test_stat_2,
  adjusted_p_value = psi_p_value_2
)
write.csv(psi_tests, output_file("psi_pairwise_tests.csv"), row.names = FALSE)

p <- ggplot(data.frame(index = seq_len(n_comparisons), statistic = psi_test_stat),
            aes(x = index, y = statistic)) +
  geom_point() + geom_hline(yintercept = c(-1.96, 1.96), color = "green2") + theme_bw()
ggsave(output_file("psi_test_statistics.png"), p, width = 8, height = 5, dpi = 300)

p <- ggplot(data.frame(index = seq_len(n_comparisons), p_value = psi_p_value_2),
            aes(x = index, y = p_value)) +
  geom_point() + geom_hline(yintercept = 0.05, color = "green") + theme_bw()
ggsave(output_file("psi_adjusted_p_values.png"), p, width = 8, height = 5, dpi = 300)

psi_pair_mat <- matrix(1, nrow = n_psi, ncol = n_psi,
                       dimnames = list(psi_index, psi_index))
psi_pair_mat[t(psi_combn)] <- psi_p_value_2
psi_pair_mat[t(psi_combn[2:1, , drop = FALSE])] <- psi_p_value_2
pheatmap(
  psi_pair_mat,
  color = colorRampPalette(c("white", "black"))(5),
  breaks = seq(0, max(0.05, max(psi_pair_mat)), length.out = 6),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  filename = output_file("psi_pairwise_pvalue_heatmap.png")
)

psi_posterior_comparisons <- do.call(rbind, lapply(
  seq_len(n_comparisons), function(i) {
    a <- psi_combn[1L, i]
    b <- psi_combn[2L, i]
    delta <- psi_sample[a, ] - psi_sample[b, ]
    ci <- quantile(delta, c(0.025, 0.975), names = FALSE)
    data.frame(
      parameter_1 = psi_index[a], parameter_2 = psi_index[b],
      mean_difference = mean(delta), posterior_sd = sd(delta),
      lower_95 = ci[1L], upper_95 = ci[2L],
      probability_positive = mean(delta > 0)
    )
  }
))
write.csv(psi_posterior_comparisons,
          output_file("psi_pairwise_posterior_comparisons.csv"), row.names = FALSE)

psi_decrease_order <- order(psi_post, decreasing = TRUE)
adjacent_diff_dist <- matrix(NA_real_, nrow = n_psi - 1L, ncol = n_record)
for (i in seq_len(n_psi - 1L)) {
  adjacent_diff_dist[i, ] <-
    psi_sample[psi_decrease_order[i], ] - psi_sample[psi_decrease_order[i + 1L], ]
}
adjacent_diff_mean <- rowMeans(adjacent_diff_dist)
Bayes_CI_95 <- t(apply(adjacent_diff_dist, 1, quantile, probs = c(0.025, 0.975)))
Bayes_CI_99 <- t(apply(adjacent_diff_dist, 1, quantile, probs = c(0.005, 0.995)))
ci_df <- data.frame(
  rank = factor(seq_len(n_psi - 1L)),
  parameter_1 = psi_index[head(psi_decrease_order, -1L)],
  parameter_2 = psi_index[tail(psi_decrease_order, -1L)],
  mean = adjacent_diff_mean,
  lower_95 = Bayes_CI_95[, 1], upper_95 = Bayes_CI_95[, 2],
  lower_99 = Bayes_CI_99[, 1], upper_99 = Bayes_CI_99[, 2]
)
write.csv(ci_df, output_file("psi_adjacent_credible_intervals.csv"), row.names = FALSE)

for (level in c("95", "99")) {
  p <- ggplot(ci_df, aes(x = rank, y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = .data[[paste0("lower_", level)]],
                      ymax = .data[[paste0("upper_", level)]])) +
    labs(title = paste0(level, "%"), x = "", y = "adjacent difference") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
  ggsave(
    output_file(paste0("psi_adjacent_CI_", level, ".png")),
    p, width = 8, height = 5, dpi = 300
  )
}

p <- ggplot(ci_df, aes(x = rank, y = mean, group = 1)) +
  geom_point() + geom_line() + geom_hline(yintercept = 0, color = "green") +
  labs(x = "", y = "adjacent difference") + theme_bw()
ggsave(output_file("psi_adjacent_mean_differences.png"), p, width = 8, height = 5, dpi = 300)

# Fisher's LSD analysis
psi_group <- factor(rep(seq_len(n_psi), each = n_record))
psi_values <- as.numeric(t(psi_sample))
m1 <- aov(psi_values ~ psi_group)
if (is.finite(deviance(m1)) && deviance(m1) > 0 && m1$df.residual > 0) {
  psi_LSD <- LSD.test(
    psi_values, psi_group, m1$df.residual,
    deviance(m1) / m1$df.residual,
    alpha = 0.05, p.adj = "BH"
  )
  capture.output(psi_LSD, file = output_file("psi_LSD_results.txt"))
  if (!is.null(psi_LSD$groups)) {
    write.csv(psi_LSD$groups, output_file("psi_LSD_groups.csv"))
  }
  if (!is.null(psi_LSD$comparison)) {
    write.csv(psi_LSD$comparison, output_file("psi_LSD_comparisons.csv"))
  }
}

# -----------------------------------------------------------------------------
# 7. Independent-mixture plug-in BIC for K = 6
# -----------------------------------------------------------------------------

BIC_mixture <- function(data_obs, .G, .K, .B, .n, .mu, .gamma, .sigma_sq, .S) {
  log_likelihood <- 0
  log_terms <- matrix(NA_real_, nrow = .n, ncol = .K)
  proportions_by_layer <- t(apply(.S, 1, function(labels) {
    tab <- table(factor(labels, levels = seq_len(.K)))
    as.numeric(tab) / sum(tab)
  }))

  for (b in seq_len(.B)) {
    for (k in seq_len(.K)) {
      mean_bk <- .mu[, k] + .gamma[, b]
      sd_bk <- sqrt(.sigma_sq[, b])
      log_density <- matrix(
        dnorm(as.numeric(data_obs[[b]]), mean = mean_bk, sd = sd_bk, log = TRUE),
        nrow = .G, ncol = .n
      )
      log_terms[, k] <- log(proportions_by_layer[b, k]) + colSums(log_density)
    }
    row_max <- apply(log_terms, 1, max)
    log_likelihood <- log_likelihood +
      sum(row_max) + sum(log(rowSums(exp(log_terms - row_max))))
  }

  # Component means, identifiable layer shifts, layer variances, and
  # layer-specific mixture proportions. The likelihood above does not use psi.
  # Posterior means and empirical weights give a plug-in criterion, not an MLE
  # BIC for the full spatial model.
  parameter_count <- .K * .G + (2 * .B - 1) * .G + .B * (.K - 1)
  -2 * log_likelihood + parameter_count * log(.B * .n)
}

BIC_K6 <- BIC_mixture(
  Y_obs_list, G, K, B, n, rocktype_post, gamma_post, sigma_sq_post, S_post
)
write.csv(
  data.frame(K = K, criterion = "independent_mixture_plugin_BIC", BIC = BIC_K6),
  output_file("BIC_K6.csv"), row.names = FALSE
)

message("Analysis complete. Every generated file was written to: ", output_dir)
