# ==============================================================================
# Manuscript:   Cultivation facilitates global naturalization
# Description:  Bayesian SEM linking native range size, cultivation extent, and
#               naturalization. One hurdle model per group provides three views:
#               incidence, extent, and overall naturalization success.
# Outputs:      Four fitted models, three 2 x 2 figures, and Table S2.
# ==============================================================================

# 1. Environment Setup ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  brms,        # Bayesian multivariate models
  posterior,   # Posterior draws
  tidyverse,   # Data manipulation and visualization
  tidygraph,   # Graph objects
  ggraph,      # SEM path diagrams
  patchwork,   # Combine panels
  ragg         # Export PNG figures
)

rm(list = ls())
cat("\014")

# --- 1.1 Paths and Analysis Settings ---
# Change project_dir if the data and results are stored elsewhere.
project_dir <- "your_work_path"
data_file <- file.path(
  project_dir, "20250530.native_plant.matched.cultivated_plant_DO11.csv"
)
result_dir <- file.path(project_dir, "results")
figure_dir <- file.path(result_dir, "figures")
model_file <- file.path(result_dir, "models_naturalization_success.RData")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

n_china <- 10L             # Possible cultivation regions within China
n_foreign <- 359L          # Possible cultivation regions outside China
seed_use <- 2024
total_cores <- parallel::detectCores(logical = TRUE)
group_labels <- c("All species", "Annual herbs", "Perennial herbs", "Woody plants")

# 2. Data Loading & Preprocessing ----------------------------------------------

# --- 2.1 Select Variables and Remove Missing Life Forms ---
native_flora <- read_csv(data_file, show_col_types = FALSE) %>%
  filter(!is.na(life.form.integrated)) %>%
  transmute(
    life.form = life.form.integrated,
    native.range = as.integer(native.global.tdwg3),
    cul.China = as.integer(planting.China.tdwg3),
    cul.foreign = as.integer(WorldCuP.n.tdwg3),
    nat.incidence = as.integer(nat.incidence),
    nat.extent = as.integer(nat.extent),
    trials.China = n_china,
    trials.foreign = n_foreign
  )

# --- 2.2 Standardize Predictors Within Each Group ---
# Apply log(x + 1) and z-score scaling separately for each life-form group.
scale_predictors <- function(data) {
  data %>%
    mutate(
      native.range.scaled = as.numeric(scale(log1p(native.range))),
      cul.China.scaled = as.numeric(scale(log1p(cul.China))),
      cul.foreign.scaled = as.numeric(scale(log1p(cul.foreign)))
    )
}

data_groups <- list(
  all = scale_predictors(native_flora),
  annual = scale_predictors(filter(native_flora, life.form == "annual herb")),
  perennial = scale_predictors(filter(native_flora, life.form == "perennial herb")),
  woody = scale_predictors(filter(native_flora, life.form == "woody"))
)

# 3. Statistical Modeling ------------------------------------------------------

# --- 3.1 Define the Bayesian SEM ---
formula_china <- bf(
  cul.China | trials(trials.China) ~ native.range.scaled,
  family = beta_binomial(link = "logit")
)

formula_foreign <- bf(
  cul.foreign | trials(trials.foreign) ~
    cul.China.scaled + native.range.scaled,
  family = beta_binomial(link = "logit")
)

formula_naturalization <- bf(
  # Positive-count component: zero-truncated negative binomial
  nat.extent ~
    cul.foreign.scaled + cul.China.scaled + native.range.scaled,
  # Hurdle component: probability of zero naturalization extent
  hu ~
    cul.foreign.scaled + cul.China.scaled + native.range.scaled,
  family = hurdle_negbinomial(link = "log", link_hu = "logit")
)

model_formula <- formula_china + formula_foreign +
  formula_naturalization + set_rescor(FALSE)

# --- 3.2 Fit the Same Model to Each Group ---
fit_success_model <- function(data) {
  brm(
    model_formula,
    data = data,
    chains = 4,
    cores = total_cores / 4,
    iter = 2000,
    seed = seed_use
  )
}

# Keep the original object names for compatibility with existing RData files.
model.success00 <- fit_success_model(data_groups$all)
model.success01 <- fit_success_model(data_groups$annual)
model.success02 <- fit_success_model(data_groups$perennial)
model.success03 <- fit_success_model(data_groups$woody)

save(
  model.success00, model.success01, model.success02, model.success03,
  file = model_file
)

# For plotting saved models only: run Section 1, load(model_file), then Sections 4-7.

# 4. Posterior Effects for Path Diagrams ---------------------------------------
models <- list(
  all = model.success00,
  annual = model.success01,
  perennial = model.success02,
  woody = model.success03
)

for (model in models) print(summary(model))

outcome_labels <- c(
  incidence = "Naturalization incidence",
  extent = "Naturalization extent",
  success = "Naturalization success"
)

outcome_scales <- c(
  incidence = "logit_P_positive",
  extent = "log_mu_positive_component",
  success = "mean_d_logEY_dx"
)

# --- 4.1 Integrate the Hurdle and Count Components for Overall Success ---
# E[Y|X] = (1 - hu) * mu / (1 - P_NB(0)).
# Each outcome edge is the mean local slope d log(E[Y|X]) / d x_j.
# These are direct associations, not total effects through cultivation mediators.
# Repeated predictor profiles are frequency-weighted; draws use batches of 128.
overall_success_draws <- function(model) {
  predictors <- c("native.range.scaled", "cul.China.scaled", "cul.foreign.scaled")
  x <- model$data[, predictors, drop = FALSE]
  if (anyNA(x) || any(!is.finite(as.matrix(x)))) {
    stop("Non-finite success predictors.")
  }

  profiles <- stats::aggregate(rep(1, nrow(x)), by = x, FUN = length)
  weights <- profiles$x / nrow(x)
  design <- t(as.matrix(profiles[, predictors, drop = FALSE]))

  coefficients <- brms::fixef(model, summary = FALSE)
  shape <- as.numeric(
    posterior::as_draws_matrix(model, variable = "shape_natextent")
  )
  mu_slopes <- coefficients[, paste0("natextent_", predictors), drop = FALSE]
  hu_slopes <- coefficients[, paste0("hu_natextent_", predictors), drop = FALSE]
  effects <- matrix(
    NA_real_, nrow(coefficients), length(predictors),
    dimnames = list(NULL, paste0("natextent_", predictors))
  )

  for (first in seq(1L, nrow(coefficients), by = 128L)) {
    ids <- first:min(first + 127L, nrow(coefficients))
    mu <- exp(sweep(
      mu_slopes[ids, , drop = FALSE] %*% design,
      1, coefficients[ids, "natextent_Intercept"], "+"
    ))
    hu <- plogis(sweep(
      hu_slopes[ids, , drop = FALSE] %*% design,
      1, coefficients[ids, "hu_natextent_Intercept"], "+"
    ))

    r <- matrix(shape[ids], nrow = length(ids), ncol = ncol(design))
    log_p0 <- -r * log1p(mu / r)
    # Derivative of log[mu / (1 - P_NB(0))] with respect to log(mu).
    sensitivity_mu <- 1 - (mu * r / (mu + r)) / expm1(-log_p0)
    mean_sensitivity <- as.vector(sensitivity_mu %*% weights)
    mean_hu <- as.vector(hu %*% weights)

    effects[ids, ] <-
      sweep(mu_slopes[ids, , drop = FALSE], 1, mean_sensitivity, "*") -
      sweep(hu_slopes[ids, , drop = FALSE], 1, mean_hu, "*")
  }

  if (any(!is.finite(effects))) {
    stop("Non-finite integrated success effects.")
  }
  return(effects)
}

# --- 4.2 Extract Draws for Each Naturalization Outcome ---
extract_path_draws <- function(model, outcome) {
  fixed_draws <- brms::fixef(model, summary = FALSE)

  if (outcome == "incidence") {
    fixed_draws <- fixed_draws[
      , str_detect(colnames(fixed_draws), "^(culChina_|culforeign_|hu_natextent_)"),
      drop = FALSE
    ]
    # logit P(Y > 0) = -logit(hu): reverse the hurdle coefficients.
    hu_columns <- str_starts(colnames(fixed_draws), "hu_natextent_")
    fixed_draws[, hu_columns] <- -fixed_draws[, hu_columns]
    colnames(fixed_draws) <- str_replace(
      colnames(fixed_draws), "^hu_natextent_", "natincidence_"
    )
  } else {
    # Extent uses log(mu) coefficients from the positive-count component.
    fixed_draws <- fixed_draws[
      , !str_starts(colnames(fixed_draws), "hu_"), drop = FALSE
    ]

    if (outcome == "success") {
      # Replace outcome slopes with integrated success effects before summarizing.
      overall_draws <- overall_success_draws(model)
      fixed_draws[, colnames(overall_draws)] <- overall_draws
      colnames(fixed_draws) <- str_replace(
        colnames(fixed_draws), "^natextent_", "natsuccess_"
      )
    }
  }

  return(fixed_draws)
}

# --- 4.3 Summarize Coefficients and Identify Diagram Edges ---
summarize_paths <- function(fixed_draws, outcome_label) {
  tibble(
    Parameter = colnames(fixed_draws),
    Median = apply(fixed_draws, 2, median),
    CI.low = apply(fixed_draws, 2, quantile, probs = 0.025),
    CI.high = apply(fixed_draws, 2, quantile, probs = 0.975),
    P.positive = colMeans(fixed_draws > 0),
    P.negative = colMeans(fixed_draws < 0)
  ) %>%
    mutate(P.direction = pmax(P.positive, P.negative)) %>%
    filter(!str_detect(Parameter, "Intercept$")) %>%
    mutate(
      to = case_when(
        str_starts(Parameter, "culChina_") ~ "Cultivation within China",
        str_starts(Parameter, "culforeign_") ~ "Cultivation outside China",
        str_detect(Parameter, "^(natincidence_|natextent_|natsuccess_)") ~
          outcome_label
      ),
      term = str_remove(
        Parameter, "^(culChina_|culforeign_|natincidence_|natextent_|natsuccess_)"
      ),
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
}

# 5. Visualization Setup -------------------------------------------------------

# --- 5.1 Shared Colors and Theme ---
color_palette <- c(
  Positive = "#C44E52",
  Negative = "#3B6FB6",
  Uncertain = "#9B9B9B"
)

theme_standard <- theme_graph(base_family = "sans") +
  theme(
    plot.title = element_text(
      size = 13, face = "bold", hjust = 0, margin = margin(b = 4)
    ),
    plot.margin = margin(8, 8, 8, 8)
  )

# --- 5.2 Draw One SEM Panel ---
# Layout: native range at the top; cultivation outside China on the left;
# cultivation within China on the right; naturalization outcome at the bottom.
plot_sem_paths <- function(path_summary, outcome_label, panel_title) {
  node_data <- tibble(
    name = c(
      "Native range size", "Cultivation within China",
      "Cultivation outside China", outcome_label
    ),
    x = c(0, 1, -1, 0),
    y = c(1.15, 0, 0, -1.15),
    display = c(
      "Native range size", "Cultivation\nwithin China",
      "Cultivation\noutside China",
      str_replace(outcome_label, "Naturalization ", "Naturalization\n")
    )
  )

  edge_data <- path_summary %>%
    mutate(
      label.position = case_when(
        from == "Cultivation within China" &
          to == "Cultivation outside China" ~ 0.42,
        from == "Native range size" & str_starts(to, "Naturalization") ~ 0.30,
        TRUE ~ 0.50
      )
    ) %>%
    select(from, to, label, Direction, Supported, label.position)

  graph <- tbl_graph(nodes = node_data, edges = edge_data, directed = TRUE)

  ggraph(graph, layout = "manual", x = x, y = y) +
    geom_edge_link(
      aes(
        label = label, colour = Direction,
        linetype = Supported, label_pos = label.position
      ),
      linewidth = 0.75,
      angle_calc = "along",
      label_dodge = grid::unit(3.2, "mm"),
      label_size = 4.2,
      label_colour = "black",
      start_cap = circle(10, "mm"),
      end_cap = rectangle(
        width = 36, height = 18, width_unit = "mm", height_unit = "mm"
      ),
      arrow = grid::arrow(length = grid::unit(2.6, "mm"), type = "closed"),
      show.legend = FALSE
    ) +
    geom_label(
      aes(x = x, y = y, label = display),
      data = node_data,
      size = 4.3,
      family = "sans",
      fontface = "bold",
      colour = "black",
      fill = "white",
      linewidth = 0.35,
      label.padding = grid::unit(2.2, "mm"),
      label.r = grid::unit(1.5, "mm")
    ) +
    scale_edge_colour_manual(values = color_palette) +
    scale_edge_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed")) +
    coord_cartesian(
      xlim = c(-1.55, 1.55), ylim = c(-1.55, 1.55), clip = "off"
    ) +
    theme_standard +
    labs(title = panel_title)
}

# 6. Generate Publication Figures ----------------------------------------------
# Each figure uses the same four groups: A = all species; B = annual herbs;
# C = perennial herbs; D = woody plants. Panel sample sizes come from model data.
figures <- list()
path_summaries <- list()

for (outcome in names(outcome_labels)) {
  outcome_label <- outcome_labels[[outcome]]
  panels <- list()

  for (i in seq_along(models)) {
    model <- models[[i]]
    fixed_draws <- extract_path_draws(model, outcome)
    path_summary <- summarize_paths(fixed_draws, outcome_label)

    panel_title <- sprintf(
      "%s   %s (n = %s)", LETTERS[i], group_labels[i],
      format(nrow(model$data), big.mark = ",", trim = TRUE)
    )
    panels[[i]] <- plot_sem_paths(path_summary, outcome_label, panel_title)

    # Keep the original Table S2 columns and the outcome/group row order.
    path_summaries[[paste(outcome, names(models)[i], sep = "_")]] <-
      path_summary %>%
      mutate(
        Figure = outcome,
        Group = group_labels[i],
        Scale = if_else(
          str_starts(to, "Naturalization"),
          outcome_scales[[outcome]], "cultivation_logit_coefficient"
        )
      )
  }

  figure <- wrap_plots(panels, ncol = 2, nrow = 2, byrow = TRUE) +
    plot_annotation(
      title = outcome_label,
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = NA)
      )
    )
  figures[[outcome]] <- figure

  ragg::agg_png(
    file.path(figure_dir, paste0("BSEM_hurdle_", outcome, ".png")),
    width = 260 / 25.4, height = 230 / 25.4, units = "in", res = 300
  )
  print(figure)
  dev.off()
}

# 7. Export Path Coefficients (Table S2) ----------------------------------------
# Three outcomes x four groups x six paths = 72 rows.
path_summary <- bind_rows(path_summaries)
write.csv(
  path_summary,
  file.path(figure_dir, "TableS2_BSEM_hurdle_coefficients.csv"),
  row.names = FALSE
)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
