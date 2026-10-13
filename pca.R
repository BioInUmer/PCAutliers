# ── Load libraries ─────────────────────────────────────────────────────
library(ggplot2)
library(MASS)
library(GGally)
library(ggnewscale)
library(patchwork)


###############################################################################
# 1. LOAD & PREPARE DATA
###############################################################################

setwd("~/Desktop/pca/")

# ── Read Eigenval ────────────────────────────────────────────────────────────
eig <- scan("1kg-pca.eigenval")
var <- round((eig / sum(eig)) * 100, 2)

pc_labels <- paste0("PC", seq_along(var))
names(var) <- pc_labels

# ── Read scores ──────────────────────────────────────────────────────────────
ref <- read.delim("1kg-projected.sscore")
names(ref)[names(ref) == "X.IID"] <- "IID"
ref$FID <- ref$IID

proj <- read.delim("target-projected.sscore")
names(proj)[names(proj) == "X.IID"] <- "IID"
proj$FID <- proj$IID

# ── Rename PC columns (flip PC1 sign) ───────────────────────────────────────
pc_cols <- paste0("PC", 1:10)

ref$PC1 <- -ref$PC1_AVG
ref$PC2 <- ref$PC2_AVG
for (i in 3:10) ref[[paste0("PC", i)]] <- ref[[paste0("PC", i, "_AVG")]]

proj$PC1 <- -proj$PC1_AVG
proj$PC2 <- proj$PC2_AVG
for (i in 3:10) proj[[paste0("PC", i)]] <- proj[[paste0("PC", i, "_AVG")]]

ref$group <- "Reference"
proj$group <- "Projected"

# ── Merge with 1KGP pedigree & assign super-populations ─────────────────────
ped <- read.delim("20130606_g1k.ped")
ped_subset <- ped[, c("Individual.ID", "Population")]

# ── Super-population mapping ────────────────────────────────────────────────
sp_levels <- c(
  "African (AFR)", "Admixed American (AMR)", "East Asian (EAS)",
  "European (EUR)", "South Asian (SAS)", "TARGET"
)

assign_superpop <- function(pop) {
  sp <- rep(NA_character_, length(pop))
  sp[pop %in% c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI")] <- "African (AFR)"
  sp[pop %in% c("CLM", "MXL", "PEL", "PUR")] <- "Admixed American (AMR)"
  sp[pop %in% c("CDX", "CHB", "CHS", "JPT", "KHV")] <- "East Asian (EAS)"
  sp[pop %in% c("CEU", "FIN", "GBR", "IBS", "TSI")] <- "European (EUR)"
  sp[pop %in% c("BEB", "GIH", "ITU", "PJL", "STU")] <- "South Asian (SAS)"
  sp[pop == "TARGET"] <- "TARGET"
  factor(sp, levels = sp_levels)
}

all_merged <- merge(
  rbind(
    ref[, c("FID", "IID", pc_cols, "group")],
    proj[, c("FID", "IID", pc_cols, "group")]
  ),
  ped_subset,
  by.x = "IID", by.y = "Individual.ID", all.x = TRUE
)
all_merged$Population[is.na(all_merged$Population)] <- "TARGET"
all_merged$SuperPop <- assign_superpop(all_merged$Population)

# Plot TARGET on top
all_merged <- all_merged[order(all_merged$SuperPop == "TARGET"), ]


###############################################################################
# 2. THEME, PALETTE & HELPERS
###############################################################################

# ── Colour palette ──────────────────────────────────────────────────────────
super_colors <- c(
  "African (AFR)"          = "#FFC20A",
  "Admixed American (AMR)" = "#5D3A9B",
  "East Asian (EAS)"       = "#E66100",
  "European (EUR)"         = "#0C7BDC",
  "South Asian (SAS)"      = "#000000",
  "TARGET"                 = "#D41159"
)

# ── Shared ggplot theme ─────────────────────────────────────────────────────
theme_pca <- function(...) {
  theme_classic(base_size = 13, base_family = "Helvetica") %+replace%
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "#333333", linewidth = 0.5),
      axis.ticks = element_line(color = "#333333", linewidth = 0.4),
      axis.text = element_text(color = "#444444", size = 10),
      axis.title.x = element_text(color = "#222222", size = 11, face = "bold", margin = margin(t = 12)),
      axis.title.y = element_text(color = "#222222", size = 11, face = "bold", margin = margin(r = 12), angle = 90),
      plot.title = element_text(
        face = "bold", size = 14, color = "#111111",
        margin = margin(b = 7)
      ),
      plot.subtitle = element_text(
        size = 9, color = "#666666",
        margin = margin(b = 10)
      ),
      plot.margin = margin(15, 15, 12, 12),
      legend.text = element_text(size = 9, color = "#333333"),
      panel.border = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      ...
    )
}

# Legend inside bottom-right
theme_legend_inset <- function() {
  theme(
    legend.position = c(0.97, 0.035),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      color = "#cccccc", linewidth = 0.35
    ),
    legend.margin = margin(6, 8, 6, 8),
    legend.key.size = unit(0.55, "cm"),
    legend.title = element_blank()
  )
}


###############################################################################
# 3. PARALLEL COORDINATES — PCs 1–10 by Ancestry
#    First plot: shows PC1–3 separate populations well, beyond PC3 they mix
###############################################################################

pc_para <- all_merged[all_merged$SuperPop != "TARGET", c("IID", pc_cols, "SuperPop")]

p_para <- ggparcoord(
  data         = pc_para,
  columns      = which(names(pc_para) %in% pc_cols),
  groupColumn  = "SuperPop",
  scale        = "uniminmax",
  alpha        = 0.25,
  splineFactor = 0
) +
  scale_color_manual(values = super_colors, name = NULL) +
  scale_x_discrete(labels = pc_cols) +
  labs(
    title = "Parallel Coordinates (1000 Genomes)",
    subtitle = "Each line = one sample · y rescaled [0,1] per PC",
    x = "PC", y = ""
  ) +
  theme_pca(
    panel.grid.major.x = element_line(color = "#dddddd", linewidth = 0.5),
    axis.ticks.y       = element_blank(),
    axis.text.y        = element_blank(),
    axis.text.x        = element_text(color = "#333333", size = 10, face = "bold"),
    legend.position    = "right"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 2.5, alpha = 1, linewidth = 1.5)
  ))


p_para
ggsave("plots/parallel-10pcs.png", p_para, width = 12, height = 6, dpi = 300, bg = "white")


###############################################################################
# 4. PCA SCATTER PLOTS — REFERENCE ONLY + REFERENCE WITH TARGET OVERLAY
###############################################################################

# ── Helper function for PCA scatter plots ───────────────────────────────────
make_pca_scatter <- function(df, xpc, ypc, var, title = NULL, subtitle = NULL) {
  ggplot(df, aes(x = .data[[xpc]], y = .data[[ypc]], color = SuperPop, shape = SuperPop)) +
    geom_point(
      data = subset(df, SuperPop != "TARGET"),
      size = 2, alpha = 0.8
    ) +
    geom_point(
      data = subset(df, SuperPop == "TARGET"),
      size = 1.5
    ) +
    scale_color_manual(values = super_colors, name = NULL) +
    scale_shape_manual(
      values = c(
        "African (AFR)"          = 16,
        "Admixed American (AMR)" = 16,
        "East Asian (EAS)"       = 16,
        "European (EUR)"         = 16,
        "South Asian (SAS)"      = 16,
        "TARGET"                 = 15
      ),
      name = NULL
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0(xpc, " (", var[xpc], "%)"),
      y = paste0(ypc, " (", var[ypc], "%)")
    ) +
    theme_pca() +
    theme_legend_inset() +
    guides(
      color = guide_legend(nrow = 1, override.aes = list(size = 2.5, alpha = 1)),
      shape = guide_legend(nrow = 1, override.aes = list(size = 2.5, alpha = 1))
    ) +
    theme(
      legend.key.size = unit(0.45, "cm")
    )
}

# ── Data subsets ────────────────────────────────────────────────────────────
ref_only <- subset(all_merged, group == "Reference")
ref_plus_target <- all_merged

# ── Reference-only PCA plots ────────────────────────────────────────────────
p_ref_12 <- make_pca_scatter(ref_only, "PC1", "PC2", var)
p_ref_13 <- make_pca_scatter(ref_only, "PC1", "PC3", var)

p_ref_combined <- p_ref_12 + p_ref_13 +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Ancestry PCA",
    subtitle = "1000 Genomes Reference",
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 16, color = "#111111",
        hjust = 0.5, margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9, color = "#666666",
        hjust = 0.5, margin = margin(b = 10)
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 15)
  )

p_ref_combined
ggsave("plots/pca-reference.png",
       p_ref_combined, width = 14, height = 7, dpi = 300, bg = "white")


# ── Reference + TARGET overlay PCA plots ────────────────────────────────────
p_target_12 <- make_pca_scatter(ref_plus_target, "PC1", "PC2", var)
p_target_13 <- make_pca_scatter(ref_plus_target, "PC1", "PC3", var)

p_target_combined <- p_target_12 + p_target_13 +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Ancestry PCA",
    subtitle = "TARGET samples projected on 1000 Genomes",
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 16, color = "#111111",
        hjust = 0.5, margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9, color = "#666666",
        hjust = 0.5, margin = margin(b = 10)
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 15)
  )

p_target_combined
ggsave("plots/pca-target.svg",
       p_target_combined, width = 14, height = 7, dpi = 300, bg = "white")


###############################################################################
# 5. MAHALANOBIS OUTLIER DETECTION — PC1 + PC2 ONLY
###############################################################################

# ── PCs used for Mahalanobis distance ───────────────────────────────────────
pc_maha <- c("PC1", "PC2")
chi2_thresh <- qchisq(p = 0.99, df = length(pc_maha))

# ── Data subsets ────────────────────────────────────────────────────────────
ref_eur <- subset(all_merged, group == "Reference" & SuperPop == "European (EUR)")
target  <- subset(all_merged, group == "Projected")


# ────────────────────────────────────────────────────────────────────────────
# 5a. APPROACH 1: EUR reference defines the center, covariance and threshold
#     Interpretation:
#     A TARGET sample is excluded if it falls outside the EUR reference ellipse.
# ────────────────────────────────────────────────────────────────────────────

eur_center <- colMeans(ref_eur[, pc_maha])
eur_cov    <- cov(ref_eur[, pc_maha])

ref_eur$maha_ref_eur <- mahalanobis(
  ref_eur[, pc_maha],
  center = eur_center,
  cov    = eur_cov
)

target$maha_ref_eur <- mahalanobis(
  target[, pc_maha],
  center = eur_center,
  cov    = eur_cov
)

target$is_outlier_ref_eur <- target$maha_ref_eur > chi2_thresh


# ────────────────────────────────────────────────────────────────────────────
# 5b. APPROACH 2: TARGET samples define their own ellipse
#     Interpretation:
#     This is fully target-driven but can be distorted by strong outliers.
# ────────────────────────────────────────────────────────────────────────────

target_center <- colMeans(target[, pc_maha])
target_cov    <- cov(target[, pc_maha])

target$maha_target <- mahalanobis(
  target[, pc_maha],
  center = target_center,
  cov    = target_cov
)

target$is_outlier_target <- target$maha_target > chi2_thresh

# ── Console summary ─────────────────────────────────────────────────────────
cat(sprintf(
  "\nEUR-reference approach: %d / %d TARGET samples flagged as outliers\n",
  sum(target$is_outlier_ref_eur), nrow(target)
))

cat(sprintf(
  "TARGET-self approach: %d / %d TARGET samples flagged as outliers\n",
  sum(target$is_outlier_target), nrow(target)
))


###############################################################################
# 6. CHI-SQUARE DIAGNOSTIC PLOTS
###############################################################################

# ── Helper function: χ² diagnostic plot ─────────────────────────────────────
make_chi2_plot <- function(df, distance_col, title = NULL, subtitle = NULL, max_x = NULL) {
  
  df <- df[is.finite(df[[distance_col]]), ]
  
  if (!is.null(max_x)) {
    n_excluded <- sum(df[[distance_col]] > max_x)
    df <- df[df[[distance_col]] <= max_x, ]
    if (n_excluded > 0) {
      subtitle_add <- sprintf("(Excluded %d extreme outliers > %g)", n_excluded, max_x)
      subtitle <- if (is.null(subtitle)) subtitle_add else paste(subtitle, subtitle_add, sep = "\n")
    }
  }
  
  ggplot(df, aes(x = .data[[distance_col]])) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth  = 1,
      boundary  = 0,
      fill      = "#5D3A9B",
      color     = "white",
      linewidth = 0.25,
      alpha     = 0.6
    ) +
    stat_function(
      aes(color = "\u03c7\u00b2(df = 2)"),
      fun       = dchisq,
      args      = list(df = length(pc_maha)),
      n         = 2000,
      linewidth = 0.8,
      lineend   = "round"
    ) +
    geom_vline(
      xintercept = chi2_thresh,
      color      = "#CC0000",
      linetype   = "dashed",
      linewidth  = 0.8
    ) +
    coord_cartesian(ylim = c(0, 0.55)) +
    scale_color_manual(values = c("\u03c7\u00b2(df = 2)" = "#CC0000")) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Mahalanobis distance²",
      y        = "Density"
    ) +
    theme_pca() +
    theme_legend_inset()
}

# ── EUR-reference χ² diagnostic ─────────────────────────────────────────────
chi2_ref_eur <- make_chi2_plot(
  df           = ref_eur,
  distance_col = "maha_ref_eur"
)

# ── TARGET-self χ² diagnostics ──────────────────────────────────────────────
chi2_target <- make_chi2_plot(
  df           = target,
  distance_col = "maha_target",
  max_x        = 50
)

chi2_combined <- chi2_ref_eur + chi2_target +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Mahalanobis Distance vs χ² Distribution",
    subtitle = "EUR Reference vs TARGET-self",
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 16, color = "#111111",
        hjust = 0.5, margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9, color = "#666666",
        hjust = 0.5, margin = margin(b = 10)
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 15)
  )

chi2_combined
ggsave("plots/chi2-maha-combined.svg",
       chi2_combined, width = 14, height = 7, dpi = 300, bg = "white")

###############################################################################
# 6.5 QQ DIAGNOSTIC PLOTS
###############################################################################

# ── Helper function: QQ plot against χ²(df = 2) ─────────────────────────────
make_qq_plot <- function(df, distance_col, title = NULL, subtitle = NULL, max_x = NULL) {
  
  df <- df[is.finite(df[[distance_col]]), ]
  
  if (!is.null(max_x)) {
    n_excluded <- sum(df[[distance_col]] > max_x)
    df <- df[df[[distance_col]] <= max_x, ]
    if (n_excluded > 0) {
      subtitle_add <- sprintf("(Excluded %d extreme outliers > %g)", n_excluded, max_x)
      subtitle <- if (is.null(subtitle)) subtitle_add else paste(subtitle, subtitle_add, sep = "\n")
    }
  }
  
  ggplot(df, aes(sample = .data[[distance_col]])) +
    stat_qq(
      distribution = stats::qchisq,
      dparams      = list(df = length(pc_maha)),
      color        = "#5D3A9B",
      size         = 1.5,
      alpha        = 0.6
    ) +
    geom_abline(
      intercept = 0,
      slope     = 1,
      color     = "#CC0000",
      linewidth = 0.8,
      linetype  = "dashed"
    ) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = sprintf("Theoretical Quantiles (\u03c7\u00b2 df = %d)", length(pc_maha)),
      y        = "Sample Quantiles (Mahalanobis distance²)"
    ) +
    theme_pca()
}

# ── EUR-reference QQ plot ───────────────────────────────────────────────────
qq_ref_eur <- make_qq_plot(
  df           = ref_eur,
  distance_col = "maha_ref_eur"
)

# ── TARGET-self QQ plot ─────────────────────────────────────────────────────
qq_target <- make_qq_plot(
  df           = target,
  distance_col = "maha_target",
  max_x        = 50
)

qq_combined <- qq_ref_eur + qq_target +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Q-Q Plot: Mahalanobis Distance vs. χ² Distribution",
    subtitle = "Comparison of EUR Reference and TARGET-self approaches",
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 14, color = "#111111",
        hjust = 0.5, margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9, color = "#666666",
        hjust = 0.5, margin = margin(b = 10)
      )
    )
  )

qq_combined
ggsave("plots/qq-maha-combined.svg",
       qq_combined, width = 14, height = 7, dpi = 300, bg = "white")


###############################################################################
# 7. ELLIPSE HELPERS AND OUTLIER PLOTS
###############################################################################

# ── Build Mahalanobis ellipse on PC1/PC2 plane ──────────────────────────────
make_maha_ellipse <- function(center, cov_mat, threshold) {
  cov_2d <- cov_mat[c("PC1", "PC2"), c("PC1", "PC2")]
  cov_eig <- eigen(cov_2d)
  
  axes <- sqrt(cov_eig$values * threshold)
  theta <- seq(0, 2 * pi, length.out = 360)
  
  ell_pts <- t(cov_eig$vectors %*% (axes * rbind(cos(theta), sin(theta))))
  
  data.frame(
    PC1 = ell_pts[, 1] + center["PC1"],
    PC2 = ell_pts[, 2] + center["PC2"]
  )
}

ellipse_ref_eur      <- make_maha_ellipse(eur_center, eur_cov, chi2_thresh)
ellipse_target       <- make_maha_ellipse(target_center, target_cov, chi2_thresh)


# ── Helper function for Mahalanobis PCA plots ───────────────────────────────
make_outlier_plot <- function(target_df, ellipse_df, outlier_col, title = NULL, subtitle = NULL, show_target = TRUE) {
  p <- ggplot() +
    geom_point(
      data = subset(all_merged, group == "Reference"),
      aes(x = PC1, y = PC2, color = SuperPop),
      size = 2, shape = 16, alpha = 0.75
    ) +
    geom_polygon(
      data = ellipse_df,
      aes(x = PC1, y = PC2),
      fill = "#004D40", color = "#004D40",
      alpha = 0.10, linewidth = 0.7
    )
    
  if (show_target) {
    p <- p +
      geom_point(
        data = target_df[!target_df[[outlier_col]], ],
        aes(x = PC1, y = PC2),
        color = "#40B0A6", size = 1.5, shape = 15
      ) +
      geom_point(
        data = target_df[target_df[[outlier_col]], ],
        aes(x = PC1, y = PC2),
        color = "#D41159", size = 1.8, shape = 17
      )
  }
  
  p +
    scale_color_manual(values = super_colors, name = NULL) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0("PC1 (", var["PC1"], "%)"),
      y = paste0("PC2 (", var["PC2"], "%)")
    ) +
    theme_pca() +
    theme_legend_inset()
}


# ── 7a. EUR-reference Mahalanobis plots ─────────────────────────────────────
p_out_ref_eur_no_target <- make_outlier_plot(
  target_df   = target,
  ellipse_df  = ellipse_ref_eur,
  outlier_col = "is_outlier_ref_eur",
  show_target = FALSE
)

p_out_ref_eur <- make_outlier_plot(
  target_df   = target,
  ellipse_df  = ellipse_ref_eur,
  outlier_col = "is_outlier_ref_eur"
)

p_out_eur_combined <- p_out_ref_eur_no_target + p_out_ref_eur +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Mahalanobis Outlier Detection (EUR Reference)",
    subtitle = sprintf("χ² df=%d, p=0.99 cutoff = %.2f", length(pc_maha), chi2_thresh),
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 16, color = "#111111",
        hjust = 0.5, margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9, color = "#666666",
        hjust = 0.5, margin = margin(b = 10)
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 15)
  )

p_out_eur_combined
ggsave("plots/pca-outliers-eur.svg",
       p_out_eur_combined, width = 14, height = 7, dpi = 300, bg = "white")


# ── 7b. TARGET-self Mahalanobis plots ───────────────────────────────────────
p_out_target_no_target <- make_outlier_plot(
  target_df   = target,
  ellipse_df  = ellipse_target,
  outlier_col = "is_outlier_target",
  show_target = FALSE
)

p_out_target <- make_outlier_plot(
  target_df   = target,
  ellipse_df  = ellipse_target,
  outlier_col = "is_outlier_target"
)

p_out_target_combined <- p_out_target_no_target + p_out_target +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Mahalanobis Outlier Detection (TARGET-self)",
    subtitle = sprintf("χ² df=%d, p=0.99 cutoff = %.2f", length(pc_maha), chi2_thresh),
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 16, color = "#111111",
        hjust = 0.5, margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9, color = "#666666",
        hjust = 0.5, margin = margin(b = 10)
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 15)
  )

p_out_target_combined
ggsave("plots/pca-outliers-target.svg",
       p_out_target_combined, width = 14, height = 7, dpi = 300, bg = "white")


# ── 7c. Combined Comparison Plot ────────────────────────────────────────────
p_out_ref_eur_orig <- make_outlier_plot(
  target_df   = target,
  ellipse_df  = ellipse_ref_eur,
  outlier_col = "is_outlier_ref_eur",
  title       = "EUR Reference"
)

p_out_target_orig <- make_outlier_plot(
  target_df   = target,
  ellipse_df  = ellipse_target,
  outlier_col = "is_outlier_target",
  title       = "TARGET-self"
)

p_out_combined <- p_out_ref_eur_orig + p_out_target_orig +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Mahalanobis Outlier Detection",
    subtitle = sprintf("χ² df=%d, p=0.99 cutoff = %.2f", length(pc_maha), chi2_thresh),
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 16, color = "#111111",
        hjust = 0.5, margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9, color = "#666666",
        hjust = 0.5, margin = margin(b = 10)
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 15)
  )

p_out_combined
ggsave("plots/pca-outliers-combined.svg",
       p_out_combined, width = 14, height = 7, dpi = 300, bg = "white")


###############################################################################
# 8. MAHALANOBIS LOLLIPOP PLOTS
###############################################################################

# ── Helper function: ranked Mahalanobis distance lollipop plot ──────────────
make_lollipop_plot <- function(df, distance_col, outlier_col, title = NULL, subtitle = NULL) {
  
  # Order samples by Mahalanobis distance
  df_ordered <- df[order(df[[distance_col]]), ]
  df_ordered$sample_rank <- seq_len(nrow(df_ordered))
  
  # Create clean labels for plotting
  df_ordered$outlier_status <- factor(
    ifelse(df_ordered[[outlier_col]], "Outlier", "Pass"),
    levels = c("Pass", "Outlier")
  )
  
  # Avoid drawing from zero on a log10 axis
  min_positive <- min(df_ordered[[distance_col]][df_ordered[[distance_col]] > 0], na.rm = TRUE)
  y_start <- min_positive / 2
  
  ggplot(
    df_ordered,
    aes(x = sample_rank, y = .data[[distance_col]], color = outlier_status)
  ) +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = chi2_thresh, ymax = Inf,
      fill = "#D41159", alpha = 0.04
    ) +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = y_start, ymax = chi2_thresh,
      fill = "#40B0A6", alpha = 0.04
    ) +
    geom_hline(
      yintercept = chi2_thresh,
      color = "#D41159",
      linetype = "dashed",
      linewidth = 0.9
    ) +
    geom_segment(
      aes(xend = sample_rank, y = y_start, yend = .data[[distance_col]]),
      linewidth = 0.3,
      alpha = 0.5
    ) +
    geom_point(size = 1.8, alpha = 0.85) +
    annotate(
      "label",
      x = nrow(df_ordered) * 0.02,
      y = chi2_thresh,
      label = sprintf(
        "χ² cutoff = %.2f\n(df=%d, p=0.99)",
        chi2_thresh, length(pc_maha)
      ),
      color = "#CC0000",
      fill = "white",
      hjust = 0,
      vjust = -0.4,
      size = 3
    ) +
    scale_y_log10(
      breaks = c(1, 5, chi2_thresh, 50, 100, 500, 1000, 5000),
      labels = function(x) round(x, 1)
    ) +
    scale_color_manual(
      values = c("Pass" = "#40B0A6", "Outlier" = "#D41159"),
      labels = c(
        "Pass"    = sprintf("Pass (n=%d)", sum(!df[[outlier_col]])),
        "Outlier" = sprintf("Outlier (n=%d)", sum(df[[outlier_col]]))
      ),
      name = NULL
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Samples ranked by Mahalanobis distance",
      y = "Mahalanobis distance² (log10)"
    ) +
    theme_pca(
      panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.4),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = c(0.02, 0.97),
      legend.justification = c("left", "top"),
      legend.background = element_rect(
        fill = alpha("white", 0.85),
        color = "#cccccc",
        linewidth = 0.35
      ),
      legend.margin = margin(6, 8, 6, 8)
    )
}


# ── EUR-reference lollipop plot ─────────────────────────────────────────────
lollipop_ref_eur <- make_lollipop_plot(
  df           = target,
  distance_col = "maha_ref_eur",
  outlier_col  = "is_outlier_ref_eur"
)

# ── TARGET-self lollipop plots ──────────────────────────────────────────────
lollipop_target <- make_lollipop_plot(
  df           = target,
  distance_col = "maha_target",
  outlier_col  = "is_outlier_target"
)

lollipop_combined <- lollipop_ref_eur + lollipop_target +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Mahalanobis Distance per Sample",
    subtitle = "Comparison of EUR Reference and TARGET-self approaches",
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 14, color = "#111111",
        hjust = 0.5, margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9, color = "#666666",
        hjust = 0.5, margin = margin(b = 10)
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 15)
  )

lollipop_combined
ggsave("plots/maha-lollipop-combined.svg",
       lollipop_combined, width = 14, height = 7, dpi = 300, bg = "white")


###############################################################################
# 9. EXPORT OUTLIERS
###############################################################################

# ── 9a. Samples excluded using EUR-reference Mahalanobis approach ───────────
outliers_ref_eur <- subset(target, is_outlier_ref_eur)[, c("FID", "IID")]

write.table(
  outliers_ref_eur,
  "samples-by-pca-ref-eur.txt",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)


# ── 9b. Samples excluded using TARGET-self Mahalanobis approach ─────────────
outliers_target_self <- subset(target, is_outlier_target)[, c("FID", "IID")]

write.table(
  outliers_target_self,
  "samples-by-pca-target-self.txt",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)