# qRT-PCR figure 2.5.1: log2 ddCt bar plots with paired t-tests for LYVE1, F13A1 and MRC1.
library(ggplot2)
library(rstatix)
library(ggpubr)
library(dplyr)
library(patchwork)

# ---- data: rows = donors, cols = Day0 / Day5_MDI / Day5_MDIX ----
LYVE1 <- matrix(c(
  0.00377064226915448, 0.004181,   3.271899,
  0.000123874512479221, 0.00012,   5.781357,
  0.000666,            0.44473,    9.155596,
  0.002144,            0.019506,   2.03257,
  0,                   2.778717,   9
), ncol = 3, byrow = TRUE)

F13A1 <- matrix(c(
  0.026045393, 0.013725666, 2.034171053,
  0.048097109, 0.054226956, 0.873293798,
  0.003285214, 0.13338127,  0.394889524
), ncol = 3, byrow = TRUE)

MRC1 <- matrix(c(
  0.05173857, 3.42539414, 150.340928,
  0.48310176, 17.6704883, 43.0116031,
  0.59177965, 15.1697457, 12.4872077
), ncol = 3, byrow = TRUE)

# ---- builder ----
build <- function(m, gene, color_dots = TRUE, correction = "none") {
  n <- nrow(m)
  df <- data.frame(
    donor = factor(rep(paste0("D", seq_len(n)), 3)),
    group = factor(rep(c("Day0", "Day5_MDI", "Day5_MDIX"), each = n),
                   levels = c("Day0", "Day5_MDI", "Day5_MDIX")),
    value = as.vector(m)
  )
  df$logval <- log2(pmax(df$value, min(df$value[df$value > 0]) / 2))  # floor zeros
  
  post <- df %>%
    pairwise_t_test(logval ~ group, paired = TRUE, p.adjust.method = correction) %>%
    add_significance(if (correction == "none") "p" else "p.adj")
  post$y.position <- max(df$value) * c(1.05, 1.18, 1.31)
  lab <- if (correction == "none") "p.signif" else "p.adj.signif"
  
  p <- ggplot(df, aes(group, value)) +
    stat_summary(fun = mean, geom = "bar", fill = "grey80", width = 0.6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    stat_pvalue_manual(post, label = lab, tip.length = 0.01) +
    theme_classic(base_size = 14) +
    labs(x = paste(gene, "expression"),
         y = expression(2^{-Delta*Delta*Ct}))
  
  if (color_dots)
    p + geom_jitter(aes(colour = donor), width = 0.1, size = 3) +
    scale_colour_brewer(palette = "Set1") + labs(colour = "Donor")
  else
    p + geom_jitter(width = 0.1, size = 3, colour = "grey25")
}

# ---- individual figures ----
ggsave("qpcr_LYVE1.png", build(LYVE1, "LYVE1"), width = 5.2, height = 4.6, dpi = 300)
ggsave("qpcr_F13A1.png", build(F13A1, "F13A1"), width = 5.2, height = 4.6, dpi = 300)
ggsave("qpcr_MRC1.png",  build(MRC1,  "MRC1"),  width = 5.2, height = 4.6, dpi = 300)

# ---- combined: neutral dots, no shared donor legend (donors not matched across genes) ----
combined <- build(LYVE1, "LYVE1", FALSE) |
  build(F13A1, "F13A1", FALSE) |
  build(MRC1,  "MRC1",  FALSE)
ggsave("qpcr_combined.png", combined, width = 13, height = 4.6, dpi = 300)