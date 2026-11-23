# ============================================================================
# Three complementary path diagrams from ONE fitted hurdle BSEM.
# Three separate 2 x 2 figures: incidence, extent, overall hurdle success.
# A = All species; B = Annual herbs; C = Perennial herbs; D = Woody plants.
# Cultivation outside China is on the LEFT; within China is on the RIGHT.
# All plotting panels are written out explicitly for straightforward editing.
# ============================================================================

cat("\014")
rm(list = ls())

library(brms)
library(dplyr)
library(tibble)
library(stringr)
library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)
library(ggplot2)

setwd("E:/中国本地植物归化/sem分析")

data.file <- "./20250530.native_plant.matched.cultivated_plant_DO11.csv"
result.dir <- "./2026.1104.results"
figure.dir <- file.path(result.dir, "figures")
dir.create(result.dir, recursive = TRUE, showWarnings = FALSE)

# The only model file read by this script.
load(file.path(result.dir, "models_naturalization_success.RData"))

# Incidence outcome: logit P(Y>0) = -logit hu (sign reversal is essential).
# Extent outcome: log mu in the zero-truncated negative binomial component.
# Success outcome: average local slope of log E[Y|X], combining hu and mu.
# Cultivation equations are shared across all three views of each model.


# For the hurdle negative binomial, E[Y|X] = (1-hu)*mu/(1-P_NB(0)).
# Each outcome edge is the empirical mean of d log(E[Y|X])/d x_j,
# holding the other predictors fixed. These are direct associations, not
# total SEM effects through cultivation mediators or original brms slopes.
# Identical predictor profiles are frequency-weighted, using every model row.
overall_success_draws <- function(model) {
  predictors <- c("native.range.scaled", "cul.China.scaled", "cul.foreign.scaled")
  x <- model$data[, predictors, drop = FALSE]
  if (anyNA(x) || any(!is.finite(as.matrix(x)))) stop("Non-finite success predictors.")
  profiles <- stats::aggregate(rep(1, nrow(x)), by = x, FUN = length)
  weights <- profiles$x / nrow(x)
  design <- t(as.matrix(profiles[, predictors, drop = FALSE]))
  coefficients <- brms::fixef(model, summary = FALSE)
  shape <- as.numeric(posterior::as_draws_matrix(model, variable = "shape_natextent"))
  mu.slopes <- coefficients[, paste0("natextent_", predictors), drop = FALSE]
  hu.slopes <- coefficients[, paste0("hu_natextent_", predictors), drop = FALSE]
  effects <- matrix(NA_real_, nrow(coefficients), length(predictors),
                    dimnames = list(NULL, paste0("natextent_", predictors)))
  for (first in seq(1L, nrow(coefficients), by = 128L)) {
    ids <- first:min(first + 127L, nrow(coefficients))
    mu <- exp(sweep(mu.slopes[ids, , drop = FALSE] %*% design, 1,
                    coefficients[ids, "natextent_Intercept"], "+"))
    hu <- plogis(sweep(hu.slopes[ids, , drop = FALSE] %*% design, 1,
                       coefficients[ids, "hu_natextent_Intercept"], "+"))
    r <- matrix(shape[ids], nrow = length(ids), ncol = ncol(design))
    log.p0 <- -r * log1p(mu / r)
    # d log[mu/(1-P_NB(0))] / d log(mu), evaluated stably.
    sensitivity.mu <- 1 - (mu * r / (mu + r)) / expm1(-log.p0)
    mean.sensitivity <- as.vector(sensitivity.mu %*% weights)
    mean.hu <- as.vector(hu %*% weights)
    effects[ids, ] <- sweep(mu.slopes[ids, , drop = FALSE], 1, mean.sensitivity, "*") -
      sweep(hu.slopes[ids, , drop = FALSE], 1, mean.hu, "*")
  }
  if (any(!is.finite(effects))) stop("Non-finite integrated success effects.")
  effects
}


# ============================================================================
# Incidence - All species
# ============================================================================

fixed.draws.inc.all <- fixef(model.success00, summary = FALSE)
fixed.draws.inc.all <- fixed.draws.inc.all[, str_detect(colnames(fixed.draws.inc.all), "^(culChina_|culforeign_|hu_natextent_)"), drop = FALSE]
hu.columns.inc.all <- str_starts(colnames(fixed.draws.inc.all), "hu_natextent_")
fixed.draws.inc.all[, hu.columns.inc.all] <- -fixed.draws.inc.all[, hu.columns.inc.all]
colnames(fixed.draws.inc.all) <- str_replace(colnames(fixed.draws.inc.all), "^hu_natextent_", "natincidence_")

fixed.summary.inc.all <- tibble(
  Parameter = colnames(fixed.draws.inc.all),
  Median = apply(fixed.draws.inc.all, 2, median),
  CI.low = apply(fixed.draws.inc.all, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.inc.all, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.inc.all > 0),
  P.negative = colMeans(fixed.draws.inc.all < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natincidence_") ~ "Naturalization incidence"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natincidence_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.inc.all <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization incidence"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nincidence"
  )
)

edge.data.inc.all <- fixed.summary.inc.all %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.inc.all <- tbl_graph(
  nodes = node.data.inc.all,
  edges = edge.data.inc.all,
  directed = TRUE
)

plot.inc.all <- ggraph(
  graph.inc.all,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.inc.all,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("A   All species (n = 25,129)", nrow(model.success00$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Incidence - Annual herbs
# ============================================================================

fixed.draws.inc.annual <- fixef(model.success01, summary = FALSE)
fixed.draws.inc.annual <- fixed.draws.inc.annual[, str_detect(colnames(fixed.draws.inc.annual), "^(culChina_|culforeign_|hu_natextent_)"), drop = FALSE]
hu.columns.inc.annual <- str_starts(colnames(fixed.draws.inc.annual), "hu_natextent_")
fixed.draws.inc.annual[, hu.columns.inc.annual] <- -fixed.draws.inc.annual[, hu.columns.inc.annual]
colnames(fixed.draws.inc.annual) <- str_replace(colnames(fixed.draws.inc.annual), "^hu_natextent_", "natincidence_")

fixed.summary.inc.annual <- tibble(
  Parameter = colnames(fixed.draws.inc.annual),
  Median = apply(fixed.draws.inc.annual, 2, median),
  CI.low = apply(fixed.draws.inc.annual, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.inc.annual, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.inc.annual > 0),
  P.negative = colMeans(fixed.draws.inc.annual < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natincidence_") ~ "Naturalization incidence"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natincidence_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.inc.annual <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization incidence"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nincidence"
  )
)

edge.data.inc.annual <- fixed.summary.inc.annual %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.inc.annual <- tbl_graph(
  nodes = node.data.inc.annual,
  edges = edge.data.inc.annual,
  directed = TRUE
)

plot.inc.annual <- ggraph(
  graph.inc.annual,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.inc.annual,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("B   Annual herbs (n = 2,002)", nrow(model.success01$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Incidence - Perennial herbs
# ============================================================================

fixed.draws.inc.perennial <- fixef(model.success02, summary = FALSE)
fixed.draws.inc.perennial <- fixed.draws.inc.perennial[, str_detect(colnames(fixed.draws.inc.perennial), "^(culChina_|culforeign_|hu_natextent_)"), drop = FALSE]
hu.columns.inc.perennial <- str_starts(colnames(fixed.draws.inc.perennial), "hu_natextent_")
fixed.draws.inc.perennial[, hu.columns.inc.perennial] <- -fixed.draws.inc.perennial[, hu.columns.inc.perennial]
colnames(fixed.draws.inc.perennial) <- str_replace(colnames(fixed.draws.inc.perennial), "^hu_natextent_", "natincidence_")

fixed.summary.inc.perennial <- tibble(
  Parameter = colnames(fixed.draws.inc.perennial),
  Median = apply(fixed.draws.inc.perennial, 2, median),
  CI.low = apply(fixed.draws.inc.perennial, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.inc.perennial, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.inc.perennial > 0),
  P.negative = colMeans(fixed.draws.inc.perennial < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natincidence_") ~ "Naturalization incidence"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natincidence_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.inc.perennial <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization incidence"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nincidence"
  )
)

edge.data.inc.perennial <- fixed.summary.inc.perennial %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.inc.perennial <- tbl_graph(
  nodes = node.data.inc.perennial,
  edges = edge.data.inc.perennial,
  directed = TRUE
)

plot.inc.perennial <- ggraph(
  graph.inc.perennial,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.inc.perennial,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("C   Perennial herbs (n = 13,772)", nrow(model.success02$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Incidence - Woody plants
# ============================================================================

fixed.draws.inc.woody <- fixef(model.success03, summary = FALSE)
fixed.draws.inc.woody <- fixed.draws.inc.woody[, str_detect(colnames(fixed.draws.inc.woody), "^(culChina_|culforeign_|hu_natextent_)"), drop = FALSE]
hu.columns.inc.woody <- str_starts(colnames(fixed.draws.inc.woody), "hu_natextent_")
fixed.draws.inc.woody[, hu.columns.inc.woody] <- -fixed.draws.inc.woody[, hu.columns.inc.woody]
colnames(fixed.draws.inc.woody) <- str_replace(colnames(fixed.draws.inc.woody), "^hu_natextent_", "natincidence_")

fixed.summary.inc.woody <- tibble(
  Parameter = colnames(fixed.draws.inc.woody),
  Median = apply(fixed.draws.inc.woody, 2, median),
  CI.low = apply(fixed.draws.inc.woody, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.inc.woody, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.inc.woody > 0),
  P.negative = colMeans(fixed.draws.inc.woody < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natincidence_") ~ "Naturalization incidence"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natincidence_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.inc.woody <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization incidence"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nincidence"
  )
)

edge.data.inc.woody <- fixed.summary.inc.woody %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.inc.woody <- tbl_graph(
  nodes = node.data.inc.woody,
  edges = edge.data.inc.woody,
  directed = TRUE
)

plot.inc.woody <- ggraph(
  graph.inc.woody,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.inc.woody,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("D   Woody plants (n = 9,355)", nrow(model.success03$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Combine incidence: All species / Annual herbs; Perennial herbs / Woody plants.
# ============================================================================
figure.incidence <- wrap_plots(
  plot.inc.all, plot.inc.annual, plot.inc.perennial, plot.inc.woody,
  ncol = 2, nrow = 2, byrow = TRUE
) +
  plot_annotation(
    title = "Naturalization incidence",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

output.incidence <- file.path(figure.dir, "BSEM_hurdle_incidence")

ragg::agg_png(paste0(output.incidence, ".png"), width = 260/25.4, height = 230/25.4,
              units = "in", res = 300)
print(figure.incidence)
dev.off()


# ============================================================================
# Extent - All species
# ============================================================================

fixed.draws.ext.all <- fixef(model.success00, summary = FALSE)
fixed.draws.ext.all <- fixed.draws.ext.all[, !str_starts(colnames(fixed.draws.ext.all), "hu_"), drop = FALSE]

fixed.summary.ext.all <- tibble(
  Parameter = colnames(fixed.draws.ext.all),
  Median = apply(fixed.draws.ext.all, 2, median),
  CI.low = apply(fixed.draws.ext.all, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.ext.all, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.ext.all > 0),
  P.negative = colMeans(fixed.draws.ext.all < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natextent_") ~ "Naturalization extent"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natextent_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.ext.all <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization extent"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nextent"
  )
)

edge.data.ext.all <- fixed.summary.ext.all %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.ext.all <- tbl_graph(
  nodes = node.data.ext.all,
  edges = edge.data.ext.all,
  directed = TRUE
)

plot.ext.all <- ggraph(
  graph.ext.all,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.ext.all,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("A   All species (n = 25,129)", nrow(model.success00$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Extent - Annual herbs
# ============================================================================

fixed.draws.ext.annual <- fixef(model.success01, summary = FALSE)
fixed.draws.ext.annual <- fixed.draws.ext.annual[, !str_starts(colnames(fixed.draws.ext.annual), "hu_"), drop = FALSE]

fixed.summary.ext.annual <- tibble(
  Parameter = colnames(fixed.draws.ext.annual),
  Median = apply(fixed.draws.ext.annual, 2, median),
  CI.low = apply(fixed.draws.ext.annual, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.ext.annual, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.ext.annual > 0),
  P.negative = colMeans(fixed.draws.ext.annual < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natextent_") ~ "Naturalization extent"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natextent_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.ext.annual <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization extent"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nextent"
  )
)

edge.data.ext.annual <- fixed.summary.ext.annual %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.ext.annual <- tbl_graph(
  nodes = node.data.ext.annual,
  edges = edge.data.ext.annual,
  directed = TRUE
)

plot.ext.annual <- ggraph(
  graph.ext.annual,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.ext.annual,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("B   Annual herbs (n = 2,002)", nrow(model.success01$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Extent - Perennial herbs
# ============================================================================

fixed.draws.ext.perennial <- fixef(model.success02, summary = FALSE)
fixed.draws.ext.perennial <- fixed.draws.ext.perennial[, !str_starts(colnames(fixed.draws.ext.perennial), "hu_"), drop = FALSE]

fixed.summary.ext.perennial <- tibble(
  Parameter = colnames(fixed.draws.ext.perennial),
  Median = apply(fixed.draws.ext.perennial, 2, median),
  CI.low = apply(fixed.draws.ext.perennial, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.ext.perennial, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.ext.perennial > 0),
  P.negative = colMeans(fixed.draws.ext.perennial < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natextent_") ~ "Naturalization extent"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natextent_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.ext.perennial <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization extent"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nextent"
  )
)

edge.data.ext.perennial <- fixed.summary.ext.perennial %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.ext.perennial <- tbl_graph(
  nodes = node.data.ext.perennial,
  edges = edge.data.ext.perennial,
  directed = TRUE
)

plot.ext.perennial <- ggraph(
  graph.ext.perennial,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.ext.perennial,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("C   Perennial herbs (n = 13,772)", nrow(model.success02$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Extent - Woody plants
# ============================================================================

fixed.draws.ext.woody <- fixef(model.success03, summary = FALSE)
fixed.draws.ext.woody <- fixed.draws.ext.woody[, !str_starts(colnames(fixed.draws.ext.woody), "hu_"), drop = FALSE]

fixed.summary.ext.woody <- tibble(
  Parameter = colnames(fixed.draws.ext.woody),
  Median = apply(fixed.draws.ext.woody, 2, median),
  CI.low = apply(fixed.draws.ext.woody, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.ext.woody, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.ext.woody > 0),
  P.negative = colMeans(fixed.draws.ext.woody < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natextent_") ~ "Naturalization extent"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natextent_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.ext.woody <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization extent"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nextent"
  )
)

edge.data.ext.woody <- fixed.summary.ext.woody %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.ext.woody <- tbl_graph(
  nodes = node.data.ext.woody,
  edges = edge.data.ext.woody,
  directed = TRUE
)

plot.ext.woody <- ggraph(
  graph.ext.woody,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.ext.woody,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("D   Woody plants (n = 9,355)", nrow(model.success03$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Combine extent: All species / Annual herbs; Perennial herbs / Woody plants.
# ============================================================================
figure.extent <- wrap_plots(
  plot.ext.all, plot.ext.annual, plot.ext.perennial, plot.ext.woody,
  ncol = 2, nrow = 2, byrow = TRUE
) +
  plot_annotation(
    title = "Naturalization extent",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

output.extent <- file.path(figure.dir, "BSEM_hurdle_extent")

ragg::agg_png(paste0(output.extent, ".png"), width = 260/25.4, height = 230/25.4,
              units = "in", res = 300)
print(figure.extent)
dev.off()

# ============================================================================
# Success - All species
# ============================================================================

fixed.draws.suc.all <- fixef(model.success00, summary = FALSE)
# Integrate hu and mu before summarising; retain a single success node.
overall.draws.suc.all <- overall_success_draws(model.success00)
fixed.draws.suc.all <- fixed.draws.suc.all[, !str_starts(colnames(fixed.draws.suc.all), "hu_"), drop = FALSE]
fixed.draws.suc.all[, colnames(overall.draws.suc.all)] <- overall.draws.suc.all
colnames(fixed.draws.suc.all) <- str_replace(colnames(fixed.draws.suc.all), "^natextent_", "natsuccess_")


fixed.summary.suc.all <- tibble(
  Parameter = colnames(fixed.draws.suc.all),
  Median = apply(fixed.draws.suc.all, 2, median),
  CI.low = apply(fixed.draws.suc.all, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.suc.all, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.suc.all > 0),
  P.negative = colMeans(fixed.draws.suc.all < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natsuccess_") ~ "Naturalization success"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natsuccess_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.suc.all <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization success"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nsuccess"
  )
)

edge.data.suc.all <- fixed.summary.suc.all %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.suc.all <- tbl_graph(
  nodes = node.data.suc.all,
  edges = edge.data.suc.all,
  directed = TRUE
)

plot.suc.all <- ggraph(
  graph.suc.all,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.suc.all,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("A   All species (n = 25,129)", nrow(model.success00$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Success - Annual herbs
# ============================================================================

fixed.draws.suc.annual <- fixef(model.success01, summary = FALSE)
# Integrate hu and mu before summarising; retain a single success node.
overall.draws.suc.annual <- overall_success_draws(model.success01)
fixed.draws.suc.annual <- fixed.draws.suc.annual[, !str_starts(colnames(fixed.draws.suc.annual), "hu_"), drop = FALSE]
fixed.draws.suc.annual[, colnames(overall.draws.suc.annual)] <- overall.draws.suc.annual
colnames(fixed.draws.suc.annual) <- str_replace(colnames(fixed.draws.suc.annual), "^natextent_", "natsuccess_")


fixed.summary.suc.annual <- tibble(
  Parameter = colnames(fixed.draws.suc.annual),
  Median = apply(fixed.draws.suc.annual, 2, median),
  CI.low = apply(fixed.draws.suc.annual, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.suc.annual, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.suc.annual > 0),
  P.negative = colMeans(fixed.draws.suc.annual < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natsuccess_") ~ "Naturalization success"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natsuccess_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.suc.annual <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization success"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nsuccess"
  )
)

edge.data.suc.annual <- fixed.summary.suc.annual %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.suc.annual <- tbl_graph(
  nodes = node.data.suc.annual,
  edges = edge.data.suc.annual,
  directed = TRUE
)

plot.suc.annual <- ggraph(
  graph.suc.annual,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.suc.annual,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("B   Annual herbs (n = 2,002)", nrow(model.success01$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Success - Perennial herbs
# ============================================================================

fixed.draws.suc.perennial <- fixef(model.success02, summary = FALSE)
# Integrate hu and mu before summarising; retain a single success node.
overall.draws.suc.perennial <- overall_success_draws(model.success02)
fixed.draws.suc.perennial <- fixed.draws.suc.perennial[, !str_starts(colnames(fixed.draws.suc.perennial), "hu_"), drop = FALSE]
fixed.draws.suc.perennial[, colnames(overall.draws.suc.perennial)] <- overall.draws.suc.perennial
colnames(fixed.draws.suc.perennial) <- str_replace(colnames(fixed.draws.suc.perennial), "^natextent_", "natsuccess_")


fixed.summary.suc.perennial <- tibble(
  Parameter = colnames(fixed.draws.suc.perennial),
  Median = apply(fixed.draws.suc.perennial, 2, median),
  CI.low = apply(fixed.draws.suc.perennial, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.suc.perennial, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.suc.perennial > 0),
  P.negative = colMeans(fixed.draws.suc.perennial < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natsuccess_") ~ "Naturalization success"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natsuccess_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.suc.perennial <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization success"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nsuccess"
  )
)

edge.data.suc.perennial <- fixed.summary.suc.perennial %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.suc.perennial <- tbl_graph(
  nodes = node.data.suc.perennial,
  edges = edge.data.suc.perennial,
  directed = TRUE
)

plot.suc.perennial <- ggraph(
  graph.suc.perennial,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.suc.perennial,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("C   Perennial herbs (n = 13,772)", nrow(model.success02$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Success - Woody plants
# ============================================================================

fixed.draws.suc.woody <- fixef(model.success03, summary = FALSE)
# Integrate hu and mu before summarising; retain a single success node.
overall.draws.suc.woody <- overall_success_draws(model.success03)
fixed.draws.suc.woody <- fixed.draws.suc.woody[, !str_starts(colnames(fixed.draws.suc.woody), "hu_"), drop = FALSE]
fixed.draws.suc.woody[, colnames(overall.draws.suc.woody)] <- overall.draws.suc.woody
colnames(fixed.draws.suc.woody) <- str_replace(colnames(fixed.draws.suc.woody), "^natextent_", "natsuccess_")


fixed.summary.suc.woody <- tibble(
  Parameter = colnames(fixed.draws.suc.woody),
  Median = apply(fixed.draws.suc.woody, 2, median),
  CI.low = apply(fixed.draws.suc.woody, 2, quantile, probs = 0.025),
  CI.high = apply(fixed.draws.suc.woody, 2, quantile, probs = 0.975),
  P.positive = colMeans(fixed.draws.suc.woody > 0),
  P.negative = colMeans(fixed.draws.suc.woody < 0)
) %>%
  mutate(P.direction = pmax(P.positive, P.negative)) %>%
  filter(!str_detect(Parameter, "Intercept$")) %>%
  mutate(
    to = case_when(
      str_starts(Parameter, "culChina_") ~ "Cultivation within China",
      str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
      str_starts(Parameter, "natsuccess_") ~ "Naturalization success"
    ),
    term = str_remove(Parameter, "^(culChina_|culforeign_|natsuccess_)"),
    from = case_when(
      term == "native.range.scaled" ~ "Native range size",
      term == "cul.China.scaled" ~ "Cultivation within China",
      term == "cul.foreign.scaled" ~ "Cultivation outside China"
    ),
    Supported = CI.low > 0 | CI.high < 0,
    Direction = case_when(
      Supported & Median > 0 ~ "Positive",
      Supported & Median < 0 ~ "Negative",
      TRUE ~ "Uncertain"
    ),
    Star = case_when(
      P.direction > 0.999 ~ "***",
      P.direction > 0.99 ~ "**",
      P.direction > 0.95 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(formatC(Median, digits = 2, format = "f"), Star)
  )

node.data.suc.woody <- tibble(
  name = c(
    "Native range size",
    "Cultivation within China",
    "Cultivation outside China",
    "Naturalization success"
  ),
  x = c(0, 1, -1, 0),
  y = c(1.15, 0, 0, -1.15),
  display = c(
    "Native range size",
    "Cultivation\nwithin China",
    "Cultivation\noutside China",
    "Naturalization\nsuccess"
  )
)

edge.data.suc.woody <- fixed.summary.suc.woody %>%
  mutate(
    label.position = case_when(
      from == "Cultivation within China" & to == "Cultivation outside China" ~ 0.42,
      from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
      TRUE ~ 0.50
    )
  ) %>%
  select(from, to, label, Direction, Supported, label.position)

graph.suc.woody <- tbl_graph(
  nodes = node.data.suc.woody,
  edges = edge.data.suc.woody,
  directed = TRUE
)

plot.suc.woody <- ggraph(
  graph.suc.woody,
  layout = "manual",
  x = x,
  y = y
) +
  geom_edge_link(
    aes(
      label = label,
      colour = Direction,
      linetype = Supported,
      label_pos = label.position
    ),
    linewidth = 0.75,
    angle_calc = "along",
    label_dodge = grid::unit(3.2, "mm"),
    label_size = 4.2,
    label_colour = "black",
    start_cap = circle(10, "mm"),
    end_cap = rectangle(
      width = 36,
      height = 18,
      width_unit = "mm",
      height_unit = "mm"
    ),
    arrow = grid::arrow(
      length = grid::unit(2.6, "mm"),
      type = "closed"
    ),
    show.legend = FALSE
  ) +
  geom_label(
    aes(x = x, y = y, label = display),
    data = node.data.suc.woody,
    size = 4.3,
    family = "sans",
    fontface = "bold",
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(2.2, "mm"),
    label.r = grid::unit(1.5, "mm")
  ) +
  scale_edge_colour_manual(
    values = c(
      Positive = "#C44E52",
      Negative = "#3B6FB6",
      Uncertain = "#9B9B9B"
    )
  ) +
  scale_edge_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dashed")
  ) +
  coord_cartesian(
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    clip = "off"
  ) +
  theme_graph(base_family = "sans") +
  labs(title = sprintf("D   Woody plants (n = 9,355)", nrow(model.success03$data))) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================================
# Combine success: All species / Annual herbs; Perennial herbs / Woody plants.
# ============================================================================
figure.success <- wrap_plots(
  plot.suc.all, plot.suc.annual, plot.suc.perennial, plot.suc.woody,
  ncol = 2, nrow = 2, byrow = TRUE
) +
  plot_annotation(
    title = "Naturalization success",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

output.success <- file.path(figure.dir, "BSEM_hurdle_success")



ragg::agg_png(paste0(output.success, ".png"), width = 260/25.4, height = 230/25.4,
              units = "in", res = 300)
print(figure.success)
dev.off()


# Export all 72 plotted paths, with scales identified.
path.summary <- bind_rows(
  mutate(fixed.summary.inc.all, Figure = "incidence", Group = "All species", Scale = if_else(str_starts(to, "Naturalization"), "logit_P_positive", "cultivation_logit_coefficient")),
  mutate(fixed.summary.inc.annual, Figure = "incidence", Group = "Annual herbs", Scale = if_else(str_starts(to, "Naturalization"), "logit_P_positive", "cultivation_logit_coefficient")),
  mutate(fixed.summary.inc.perennial, Figure = "incidence", Group = "Perennial herbs", Scale = if_else(str_starts(to, "Naturalization"), "logit_P_positive", "cultivation_logit_coefficient")),
  mutate(fixed.summary.inc.woody, Figure = "incidence", Group = "Woody plants", Scale = if_else(str_starts(to, "Naturalization"), "logit_P_positive", "cultivation_logit_coefficient")),
  mutate(fixed.summary.ext.all, Figure = "extent", Group = "All species", Scale = if_else(str_starts(to, "Naturalization"), "log_mu_positive_component", "cultivation_logit_coefficient")),
  mutate(fixed.summary.ext.annual, Figure = "extent", Group = "Annual herbs", Scale = if_else(str_starts(to, "Naturalization"), "log_mu_positive_component", "cultivation_logit_coefficient")),
  mutate(fixed.summary.ext.perennial, Figure = "extent", Group = "Perennial herbs", Scale = if_else(str_starts(to, "Naturalization"), "log_mu_positive_component", "cultivation_logit_coefficient")),
  mutate(fixed.summary.ext.woody, Figure = "extent", Group = "Woody plants", Scale = if_else(str_starts(to, "Naturalization"), "log_mu_positive_component", "cultivation_logit_coefficient")),
  mutate(fixed.summary.suc.all, Figure = "success", Group = "All species", Scale = if_else(str_starts(to, "Naturalization"), "mean_d_logEY_dx", "cultivation_logit_coefficient")),
  mutate(fixed.summary.suc.annual, Figure = "success", Group = "Annual herbs", Scale = if_else(str_starts(to, "Naturalization"), "mean_d_logEY_dx", "cultivation_logit_coefficient")),
  mutate(fixed.summary.suc.perennial, Figure = "success", Group = "Perennial herbs", Scale = if_else(str_starts(to, "Naturalization"), "mean_d_logEY_dx", "cultivation_logit_coefficient")),
  mutate(fixed.summary.suc.woody, Figure = "success", Group = "Woody plants", Scale = if_else(str_starts(to, "Naturalization"), "mean_d_logEY_dx", "cultivation_logit_coefficient"))
)
write.csv(path.summary, file.path(figure.dir, "TableS2_BSEM_hurdle_coefficients.csv"), row.names = FALSE)

