# ── Colour palette ──────────────────────────────────────────────────────────
super_colors <- c(
    "African (AFR)"          = "#FFC107",
    "Admixed American (AMR)" = "#5D3A9B",
    "East Asian (EAS)"       = "#E66100",
    "European (EUR)"         = "#0C7BDC",
    "South Asian (SAS)"      = "#000000",
    "TARGET"                 = "#D81B60"
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
                face = "bold", size = 17, color = "#111111",
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
