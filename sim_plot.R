#!/usr/bin/env Rscript

# Standalone version of unsupervised_K6_v2.Rmd
# Run from the directory containing alpha.txt, mu.txt, gamma.txt, etc.:
#   Rscript unsupervised_K6_v2.R -o results

parse_args <- function(args) {
  output_dir <- NULL
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript unsupervised_K6_v2.R -o OUTPUT_DIR\n\n",
        "Options:\n",
        "  -o, --output DIR   Directory for every generated file\n",
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
  data.matrix(fread(path, sep = " ", header = FALSE))
}

read_label_matrix <- function(path, expected_columns = 4600L) {
  labels <- as.matrix(read.table(path, sep = " ", header = FALSE, fill = TRUE))
  while (ncol(labels) > expected_columns && all(is.na(labels[, ncol(labels)]))) {
    labels <- labels[, -ncol(labels), drop = FALSE]
  }
  if (ncol(labels) != expected_columns) {
    stop(
      basename(path), " has ", ncol(labels), " label columns; expected ",
      expected_columns, ".", call. = FALSE
    )
  }
  labels + 1L
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
observed_file <- input_file("zhunbei_cpp_3.0.txt")
labels_file <- input_file("post_labels.txt")
depth_file <- normalizePath(file.path(input_dir, "..", "inversion_depth.dat"), mustWork = FALSE)
coordinates_file <- normalizePath(
  file.path(input_dir, "..", "inversion_coordinates.txt"), mustWork = FALSE
)

assert_files_exist(c(
  alpha_file, mu_file, gamma_file, sigma_sq_file, psi_file,
  observed_file, labels_file, depth_file, coordinates_file
))

alpha_record <- read_numeric_matrix(alpha_file)
mu_record <- read_numeric_matrix(mu_file)
gamma_record <- read_numeric_matrix(gamma_file)
sigma_sq_record <- read_numeric_matrix(sigma_sq_file)
psi_record <- read_numeric_matrix(psi_file)
Y_obs <- read_numeric_matrix(observed_file)
S_post <- read_label_matrix(labels_file)

n_iter <- 10000L
n_record <- n_iter / 2L
after_burnin <- (n_record + 1L):n_iter
G <- 3L
K <- 6L
B <- 11L
n <- 4600L

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

if (!identical(dim(Y_obs), c(G, n * B))) {
  stop(
    "zhunbei_cpp_3.0.txt must have dimensions 3 x 50,600; found ",
    paste(dim(Y_obs), collapse = " x "), ".", call. = FALSE
  )
}
if (!identical(dim(S_post), c(B, n))) {
  stop(
    "post_labels.txt must have dimensions 11 x 4,600; found ",
    paste(dim(S_post), collapse = " x "), ".", call. = FALSE
  )
}

rownames(Y_obs) <- c("log(resistivity)", "density", "velocity")
colnames(Y_obs) <- paste0(rep(seq_len(B), each = n), "_", rep(seq_len(n), B))

time_to_depth <- read.table(depth_file)
colnames(time_to_depth) <- c("time", "depth")
inversion_depth <- time_to_depth[time_to_depth$time %in% seq(3000, 3500, 50), ]
inversion_depth$depth <- round(inversion_depth$depth / 1000, digits = 2)
if (nrow(inversion_depth) != B) {
  stop("Expected 11 selected inversion depths; found ", nrow(inversion_depth), ".", call. = FALSE)
}

inversion_coordinates <- read.table(coordinates_file, sep = " ", header = FALSE)
colnames(inversion_coordinates) <- c("x", "y")
if (nrow(inversion_coordinates) != n) {
  stop("Expected 4,600 inversion coordinates; found ", nrow(inversion_coordinates), ".", call. = FALSE)
}

well_coordinates <- cbind(
  c(15570383.24, 15573693.94, 15553318.27),
  c(5079666.25, 5078575.90, 5079825.38)
) / 1000

Y_melt <- reshape2::melt(Y_obs)
colnames(Y_melt) <- c("property", "sample", "value")
Y_melt$depth <- rep(inversion_depth$depth, each = n * G)

# -----------------------------------------------------------------------------
# 2. Posterior summaries and summary CSV files
# -----------------------------------------------------------------------------

alpha_post <- rowMeans(alpha_record[, after_burnin, drop = FALSE])
mu_post <- t(matrix(
  rowMeans(mu_record[, after_burnin, drop = FALSE]), nrow = K, ncol = G
))
gamma_post <- t(matrix(
  rowMeans(gamma_record[, after_burnin, drop = FALSE]), nrow = B, ncol = G
))
sigma_sq_post <- t(matrix(
  rowMeans(sigma_sq_record[, after_burnin, drop = FALSE]), nrow = B, ncol = G
))
subtype_post <- sweep(mu_post, 1, alpha_post, FUN = "+")

names(alpha_post) <- paste0("alpha_", seq_len(G))
colnames(mu_post) <- paste0("mu_", seq_len(K))
colnames(subtype_post) <- paste0("subtype", seq_len(K))
rownames(subtype_post) <- c("R", "D", "V")
rownames(gamma_post) <- c("R", "D", "V")
colnames(gamma_post) <- paste0("gamma", seq_len(B))
rownames(sigma_sq_post) <- c("R", "D", "V")
colnames(sigma_sq_post) <- paste0("sigma_sq", seq_len(B))

write.csv(alpha_post, output_file("alpha_K6.csv"))
write.csv(mu_post, output_file("mu_K6.csv"))
write.csv(subtype_post, output_file("subtype_K6.csv"))
write.csv(t(gamma_post), output_file("gamma_K6.csv"))
write.csv(t(sigma_sq_post), output_file("sigma_sq_K6.csv"))

# -----------------------------------------------------------------------------
# 3. Parameter heatmaps
# -----------------------------------------------------------------------------

subtype_est <- mu_post
rownames(subtype_est) <- c("resistivity", "density", "velocity")
subtype_melt <- reshape2::melt(subtype_est)
colnames(subtype_melt) <- c("property", "subtype", "value")
subtype_melt$subtype <- factor(subtype_melt$subtype)

p <- ggplot(subtype_melt, aes(x = subtype, y = property, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3))) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-1.5, 1.5)
  ) +
  scale_x_discrete(position = "top") +
  labs(title = "Estimated subtype effects, K=6", x = "subtype", y = "property") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
ggsave(output_file("estimated_subtype_effects.png"), p, width = 7, height = 4, dpi = 300)

gamma_plot_matrix <- gamma_post
rownames(gamma_plot_matrix) <- c("resistivity", "density", "velocity")
colnames(gamma_plot_matrix) <- format(inversion_depth$depth, nsmall = 2)
gamma_melt <- reshape2::melt(gamma_plot_matrix)
colnames(gamma_melt) <- c("property", "depth", "value")
gamma_melt$depth <- factor(gamma_melt$depth, levels = colnames(gamma_plot_matrix))

p <- ggplot(gamma_melt, aes(x = depth, y = property, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3))) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-0.35, 0.35)
  ) +
  scale_x_discrete(position = "top") +
  labs(title = "Estimated depth effects", x = "depth (km)", y = "property") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
ggsave(output_file("estimated_depth_effects.png"), p, width = 10, height = 4, dpi = 300)

sigma_melt <- reshape2::melt(sigma_sq_post)
colnames(sigma_melt) <- c("property", "depth", "value")
p <- ggplot(sigma_melt, aes(x = depth, y = property, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  scale_x_discrete(position = "top") +
  labs(title = "Estimated variance effects", x = "depth", y = "property") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
ggsave(output_file("estimated_variance_effects.png"), p, width = 9, height = 4, dpi = 300)

# Publication-style parameter tables from the notebook
subtype_app_df <- reshape2::melt(subtype_post)
colnames(subtype_app_df) <- c("property", "lithology", "value")
subtype_app_df$lithology <- factor(subtype_app_df$lithology)
p <- ggplot(subtype_app_df, aes(x = lithology, y = property, fill = value)) +
  geom_tile() +
  geom_text(aes(label = format(round(value, 2), nsmall = 2)), size = 7, color = "green3") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "", y = "") +
  theme(
    axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 22),
    legend.text = element_text(size = 20), legend.title = element_text(size = 22)
  )
ggsave(output_file("subtype_app.png"), p, dpi = 100, width = 6, height = 3, units = "in")

gamma_app_df <- gamma_melt
p <- ggplot(gamma_app_df, aes(x = depth, y = property, fill = value)) +
  geom_tile() +
  geom_text(aes(label = format(round(value, 2), nsmall = 2)), size = 6, color = "green3") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(
    low = "skyblue2", mid = "white", high = "#E3095C",
    midpoint = 0, limits = c(-0.35, 0.35)
  ) +
  labs(x = "", y = "") +
  theme(
    axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 20),
    legend.text = element_text(size = 20), legend.title = element_text(size = 22)
  )
ggsave(output_file("gamma_app.png"), p, dpi = 100, width = 9, height = 3, units = "in")

# heatmap3 versions
subtype_heatmap <- subtype_post[3:1, , drop = FALSE]
colnames(subtype_heatmap) <- paste0("type ", seq_len(K))
rownames(subtype_heatmap) <- c("V", "D", "R")
png(output_file("est_subtype.png"), width = 600, height = 500)
heatmap3(
  subtype_heatmap, Colv = NA, Rowv = NA, scale = "none",
  ColSideColors = c("#E8585D", "cornflowerblue", "#50BF50", "orange", "pink", "yellow"),
  col = colorRampPalette(c("white", "black"))(20),
  ColSideLabs = FALSE, labRow = NA, labCol = NA
)
dev.off()

depth_heatmap <- gamma_plot_matrix[3:1, , drop = FALSE]
colnames(depth_heatmap) <- paste0("L", seq_len(B))
rownames(depth_heatmap) <- c("V", "D", "R")
gradient_colors <- colorRampPalette(c("#00dbde", "#fc00ff"))(B)
png(output_file("est_depth.png"), width = 600, height = 500)
heatmap3(
  depth_heatmap, Colv = NA, Rowv = NA, scale = "none",
  ColSideColors = gradient_colors,
  col = colorRampPalette(c("blue", "white", "red"))(20),
  ColSideLabs = FALSE, labRow = NA, labCol = NA
)
dev.off()

# -----------------------------------------------------------------------------
# 4. Spatial clustering maps
# -----------------------------------------------------------------------------

color_box <- c("#E8585D", "cornflowerblue", "#50BF50", "orange", "pink", "yellow")
well_colors <- c("red3", "blue", "green")
well_size <- 2

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

  for (j in 1:3) {
    p <- p + geom_rect(
      xmin = well_coordinates[j, 1] - well_size,
      xmax = well_coordinates[j, 1] + well_size,
      ymin = well_coordinates[j, 2] - well_size,
      ymax = well_coordinates[j, 2] + well_size,
      color = well_colors[j], fill = NA, linewidth = 1
    )
  }

  ggsave(
    output_file(paste0("Clustering_K6_", b, ".png")),
    p, units = "in", width = 5, height = 2.5, dpi = 150
  )
}

# -----------------------------------------------------------------------------
# 5. Observed-data plots by depth and predicted label
# -----------------------------------------------------------------------------

S_vec <- unlist(t(S_post))
Y_melt$label <- factor(rep(S_vec, each = G), levels = seq_len(K))

p <- ggplot(
  Y_melt[Y_melt$property == "log(resistivity)", ],
  aes(x = label, y = value, fill = factor(depth))
) +
  geom_boxplot() +
  labs(x = "predicted label", y = "log(resistivity)", fill = "depth (km)") +
  theme_bw()
ggsave(output_file("resistivity_by_label_and_depth.png"), p, width = 10, height = 6, dpi = 300)

gradient_colors <- colorRampPalette(c("#00dbde", "#fc00ff"))(B)
label1_resistivity <- Y_melt[
  Y_melt$property == "log(resistivity)" & Y_melt$label == "1", ]
label1_median <- aggregate(value ~ depth, data = label1_resistivity, FUN = median)
depth_levels <- rev(unique(as.character(inversion_depth$depth)))
label1_resistivity$depth <- factor(label1_resistivity$depth, levels = depth_levels)
label1_median$depth <- factor(label1_median$depth, levels = depth_levels)

p <- ggplot() +
  geom_boxplot(
    data = label1_resistivity,
    aes(x = depth, y = value, fill = depth), show.legend = FALSE
  ) +
  geom_line(
    data = label1_median,
    aes(x = depth, y = value, group = 1),
    color = "#FFD700", linewidth = 1.2, linetype = "longdash"
  ) +
  scale_fill_manual(values = gradient_colors) +
  scale_y_continuous(position = "right") +
  coord_flip() +
  labs(x = "depth (km)", y = "R") +
  theme_bw() +
  theme(
    panel.grid = element_blank(), panel.border = element_blank(),
    axis.line = element_line(color = "black"), axis.title = element_text(size = 20),
    axis.text = element_text(size = 14)
  )
ggsave(output_file("resistivity_label1_by_depth.png"), p, width = 7, height = 6, dpi = 300)

# Scatter plots that can be produced without the undefined S_post_3 object
scatter_labels <- factor(unlist(t(S_post)), levels = seq_len(K))
p <- ggplot(data.frame(R = Y_obs[1, ], D = Y_obs[2, ], label = scatter_labels),
            aes(x = R, y = D, color = label)) +
  geom_point(size = 0.4) +
  labs(color = "label", y = "density", x = "resistivity", title = "I1")
ggsave(output_file("resistivity_vs_density.png"), p, width = 7, height = 5, dpi = 300)

p <- ggplot(data.frame(R = Y_obs[1, ], V = Y_obs[3, ], label = scatter_labels),
            aes(x = R, y = V, color = label)) +
  geom_point(size = 0.4) +
  labs(color = "label", y = "velocity", x = "resistivity", title = "I1")
ggsave(output_file("resistivity_vs_velocity.png"), p, width = 7, height = 5, dpi = 300)

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
acceptance_rate <- rowMeans(psi_record[, -1, drop = FALSE] != psi_record[, -ncol(psi_record), drop = FALSE])
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
      (2 * n_record - 1)
  )
  psi_test_stat[i] <- mean_diff[i] / (sd_pool[i] * sqrt(2 / n_record))
  psi_test_stat_2[i] <- mean_diff[i] / sqrt(
    var(psi_sample[idx1, ]) + var(psi_sample[idx2, ])
  )
}

psi_p_value <- 2 * pt(-abs(psi_test_stat), df = 2 * n_record - 1)
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

psi_pair_mat <- matrix(0, nrow = n_psi, ncol = n_psi, dimnames = list(psi_index, psi_index))
psi_pair_mat[upper.tri(psi_pair_mat)] <- psi_p_value_2
psi_pair_mat <- psi_pair_mat + t(psi_pair_mat)
pheatmap(
  psi_pair_mat,
  color = colorRampPalette(c("white", "black"))(5),
  breaks = seq(0, max(0.05, max(psi_pair_mat)), length.out = 6),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  filename = output_file("psi_pairwise_pvalue_heatmap.png")
)

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
  geom_point() + geom_line() + geom_hline(yintercept = 0.05, color = "green") +
  labs(x = "", y = "adjacent difference") + theme_bw()
ggsave(output_file("psi_adjacent_mean_differences.png"), p, width = 8, height = 5, dpi = 300)

# Fisher's LSD analysis
psi_group <- factor(rep(seq_len(n_psi), each = n_record))
psi_values <- as.numeric(t(psi_sample))
m1 <- aov(psi_values ~ psi_group)
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

# -----------------------------------------------------------------------------
# 7. BIC for K = 6
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
      log_terms[, k] <- log(proportions_by_layer[b, k]) +
        colSums(dnorm(data_obs[, , b], mean = mean_bk, sd = sd_bk, log = TRUE))
    }
    row_max <- apply(log_terms, 1, max)
    log_likelihood <- log_likelihood +
      sum(row_max) + sum(log(rowSums(exp(log_terms - row_max))))
  }

  parameter_count <- .K * .G + (2 * .B - 1) * .G + .K * (.K - 1) / 2
  -2 * log_likelihood + parameter_count * log(.G * .B * .n)
}

Y_array <- array(as.numeric(Y_obs), dim = c(G, n, B))
BIC_K6 <- BIC_mixture(
  Y_array, G, K, B, n, subtype_post, gamma_post, sigma_sq_post, S_post
)
write.csv(data.frame(K = K, BIC = BIC_K6), output_file("BIC_K6.csv"), row.names = FALSE)

message("Analysis complete. Every generated file was written to: ", output_dir)
