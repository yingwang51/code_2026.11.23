# ============================================================================
# Four-path decomposition of naturalization incidence, extent and success
# from the same integrated success model (g-computation, tidyverse style)
# ============================================================================
#
# PURPOSE
#   For each species group (All, Annual, Perennial, Woody), this script splits
#   the effect of native range size (X) on naturalization (Y) into a direct
#   path and three indirect paths through cultivation in China (C) and
#   cultivation outside China (F). The fitted brms model is used as a
#   generative model: the mediators are simulated from their fitted
#   distributions, and the outcome is predicted under a set of
#   counterfactual scenarios (parametric g-computation).
#
# MODEL STRUCTURE (one multivariate brms model per species group)
#   C: cul.China   ~ beta-binomial(trials.China)
#        logit(p_C) = c0 + c1 * X
#   F: cul.foreign ~ beta-binomial(trials.foreign)
#        logit(p_F) = f0 + f1 * X + f2 * C
#   Y: nat.extent  ~ hurdle negative binomial
#        log(mu)    = a0 + a1 * X + a2 * C + a3 * F   (count part)
#        logit(hu)  = h0 + h1 * X + h2 * C + h3 * F   (zero part)
#   X, C and F enter the model as z-scores of log1p(count), computed within
#   each group. Predictions use fixed effects only; group-level (random)
#   effects, if present in the model, are not used.
#   X is the only exogenous variable, so no averaging over observed covariates
#   is needed: the data are used only to obtain the percentiles of X and the
#   standardization constants.
#
# CONTRAST
#   X is raised from its 25th to its 75th percentile, computed on the observed
#   (untransformed) native range values of the group. C and F are not fixed:
#   they are simulated from the fitted beta-binomial distributions given the
#   scenario values of their parents (see FOUR PATHS).
#
# THREE ESTIMANDS (hurdle negative binomial)
#   mu is the mean of the untruncated NB distribution (not the mean of Y);
#   P_NB(0) = (1 + mu / shape)^(-shape) is its probability of a zero.
#   incidence = 1 - hu
#               Pr(Y > 0): probability of naturalizing in at least one region
#   extent    = mu / [1 - P_NB(0)]
#               E[Y | Y > 0]: standardized (zero-truncated) mean of the
#               positive count component, i.e. the expected number of
#               naturalized regions given that the species has naturalized
#   success   = (1 - hu) * mu / [1 - P_NB(0)] = incidence * extent
#               E[Y]: unconditional expected number of naturalized regions
#   Each estimand is computed within every mediator simulation and then
#   averaged, first over simulations (columns) and then over posterior draws
#   (rows). Hence mean(extent) is NOT mean(success) / mean(incidence): the
#   ratio of means would weight each simulation by its incidence, whereas
#   averaging extent directly gives all simulations equal weight.
#
# FOUR PATHS (switched on in the order Direct -> China -> Foreign -> Both)
#   Notation: x0 = Q25 and x1 = Q75 of X; C(x) is C simulated at X = x;
#   F(x, c) is F simulated at X = x and C = c.
#
#   Scenario   X in Y   C in Y   F in Y
#   Baseline   x0       C(x0)    F(x0, C(x0))
#   Direct     x1       C(x0)    F(x0, C(x0))
#   China      x1       C(x1)    F(x0, C(x0))
#   Foreign    x1       C(x1)    F(x1, C(x0))
#   Both       x1       C(x1)    F(x1, C(x1))
#
#   Direct      = Direct  - Baseline   X -> Y
#   Via_China   = China   - Direct     X -> C -> Y
#   Via_Foreign = Foreign - China      X -> F -> Y
#   Via_Both    = Both    - Foreign    X -> C -> F -> Y
#   Total       = Both    - Baseline   (= sum of the four paths)
#   Because the link functions are non-linear, the size of each path depends
#   on the order in which the paths are switched on; only the total effect is
#   independent of this order.
#
# STRUCTURAL ASSUMPTION (cross-scenario coupling)
#   Each mediator is simulated by inverse-CDF sampling from one fixed set of
#   uniform random numbers that is re-used in every scenario. A simulation
#   (column) therefore represents the same latent "unit" in all scenarios:
#   its C values at x0 and x1 are perfectly rank-correlated, as are its F
#   values across the F scenarios.
#   This is an additional structural assumption, not merely a device to
#   reduce Monte Carlo error. In the China and Foreign scenarios, Y receives
#   C(x1) directly while F was generated from C(x0). The joint distribution
#   of C(x0) and C(x1) cannot be learned from the data; the coupling fixes
#   it, and therefore affects the point estimates of these two scenarios.
#   (Coupling of F only reduces the Monte Carlo error of the path effects,
#   because each scenario uses a single F.)
#   The coupling does not by itself guarantee causal identification, which
#   also requires the usual no-unmeasured-confounding assumptions.
#
# OTHER NOTES
#   - Models are not refitted; only their posterior draws are re-used.
#   - Only posterior means are reported; no credible intervals.
#   - No custom functions are defined in the g-computation part (Sections 1-8).
#
# REQUIREMENTS
#   R >= 4.0 (named error messages in stopifnot()); purrr >= 1.0
#   (list_rbind()). The pipe %>% comes from magrittr and is loaded with
#   tidyverse.
# ============================================================================

library(brms)       # fixef(), posterior_epred()
library(posterior)  # as_draws_matrix()
library(tidyverse)  # dplyr, readr, tibble, purrr, tidyr and the pipe %>%

# ============================================================================
# 1. Paths and simulation settings
# ============================================================================

# Working directory. This path contains the only non-ASCII characters in the
# script, so save the file as UTF-8. Before sharing the code, replace setwd()
# with a relative path or an RStudio project.
setwd("E:/中国本地植物归化/sem分析")
setwd("")

model_path <- "./2026.1104.results/models_naturalization_success.RData"
data_path  <- "./20250530.native_plant.matched.cultivated_plant_DO11.csv"
result_dir <- "./2026.1104.results"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Monte Carlo settings
#   Each mediator scenario is stored as an n_draws x n_sims matrix of doubles
#   (about 16 MB at the default 1000 x 2000). Several such matrices are held
#   in memory at the same time; reduce n_sims or max_draws if memory is
#   limited.
n_sims    <- 2000L  # mediator simulations per posterior draw
max_draws <- 1000L  # posterior draws used; set to Inf to use all draws
base_seed <- 2024L  # seed for group i is base_seed + i. Resetting the seed
                    # for each group makes its results reproducible on their
                    # own, provided the order of `models` is unchanged.

# Estimands and the units written to the output table. All effects are
# differences between scenario means on the response scale.
outcomes <- c("incidence", "extent", "success")
outcome_units <- c(
  incidence = "probability difference",
  extent    = "difference in standardized positive-component mean regions",
  success   = "difference in expected naturalization regions"
)

# ---- Load models: one integrated model per species group ----
# load() restores the object model.env, which holds the four fitted brmsfit
# objects. The list names below become the Group labels in all output tables.
model.env <- new.env()
load(model_path, envir = model.env)

models <- list(
  All       = model.env$model.success00,
  Annual    = model.env$model.success01,
  Perennial = model.env$model.success02,
  Woody     = model.env$model.success03
)

# ---- Read raw data ----
# fit$data stores native range only in standardized form, so the raw counts
# are needed to (i) compute percentiles on the original scale and (ii)
# recover the standardization constants. The two cultivation columns are
# read only to check that the rows match fit$data (Section 3.1).
# Species without a life-form record are dropped, as in the data used for
# fitting; Section 3.1 stops the script if the rows do not match.
raw_data <- read_csv(
  data_path,
  show_col_types = FALSE,
  col_select = c(life.form.integrated, native.global.tdwg3,
                 planting.China.tdwg3, WorldCuP.n.tdwg3)
) %>%
  filter(!is.na(life.form.integrated))

# Group subsets. Each must contain the same species, in the same order, as
# the data used to fit the corresponding model.
raw_data_by_group <- list(
  All       = raw_data,
  Annual    = filter(raw_data, life.form.integrated == "annual herb"),
  Perennial = filter(raw_data, life.form.integrated == "perennial herb"),
  Woody     = filter(raw_data, life.form.integrated == "woody")
)

# ============================================================================
# 2. Scenario definitions (shared by all groups)
#    The columns x and china take "low" or "high", i.e. the values obtained
#    with X set to Q25 (x0) or Q75 (x1). See the scenario table in the header.
# ============================================================================

# Three scenarios for foreign cultivation F = F(X, C)
#   x     : level of X in the F equation
#   china : which simulated C, C(x0) or C(x1), enters the F equation
#
#   xlow_clow   = F(x0, C(x0))   baseline
#   xhigh_clow  = F(x1, C(x0))   opens X -> F only
#   xhigh_chigh = F(x1, C(x1))   additionally opens X -> C -> F
foreign_scenarios <- tribble(
  ~name,         ~x,     ~china,
  "xlow_clow",   "low",  "low",
  "xhigh_clow",  "high", "low",
  "xhigh_chigh", "high", "high"
)

# Five scenarios for the outcome Y = Y(X, C, F)
#   x       : level of X in the Y equation
#   china   : which simulated C enters the Y equation
#   foreign : which simulated F (a name from foreign_scenarios) enters Y
# Each row opens one more path than the row above, so the difference between
# adjacent rows is one path:
#   Direct      = Direct  - Baseline
#   Via_China   = China   - Direct
#   Via_Foreign = Foreign - China
#   Via_Both    = Both    - Foreign
outcome_scenarios <- tribble(
  ~name,      ~x,     ~china, ~foreign,
  "Baseline", "low",  "low",  "xlow_clow",
  "Direct",   "high", "low",  "xlow_clow",
  "China",    "high", "high", "xlow_clow",
  "Foreign",  "high", "high", "xhigh_clow",
  "Both",     "high", "high", "xhigh_chigh"
)

# Coefficient prefixes of the two hurdle parts.
# In multivariate models, brms names fixed effects "<response>_<term>" for
# the main parameter and "<dpar>_<response>_<term>" for other distributional
# parameters, with dots removed from response names (nat.extent -> natextent,
# cul.China -> culChina).
#   mu = mean of the untruncated NB count distribution (log link)
#   hu = probability of a zero (logit link)
hurdle_prefix <- c(mu = "natextent_", hu = "hu_natextent_")

# Containers for results, keyed by "<group>_<outcome>" and combined in
# Section 8
effect_tables   <- list()  # path effects and relative contributions
scenario_tables <- list()  # posterior mean of each scenario
check_tables    <- list()  # null-effect and additivity checks

# ============================================================================
# 3. G-computation for each species group
# ============================================================================

for (i in seq_along(models)) {
  group            <- names(models)[i]
  fit              <- models[[i]]
  fit_data         <- fit$data                        # rows used in fitting
  group_raw        <- raw_data_by_group[[group]]
  native_range_raw <- group_raw$native.global.tdwg3   # raw native range (X)

  # ---- 3.1 Check that the raw subset matches the fitted data row by row ----
  # The percentiles and standardization constants below are computed from
  # group_raw, so it must contain exactly the species used in fitting, in the
  # same order. Identical cultivation counts in every row serve as a
  # fingerprint of row identity and order. stopifnot() evaluates the
  # conditions in turn, so the length check runs before the element-wise
  # comparisons.
  stopifnot(
    "Number of species differs from the fitted data" =
      length(native_range_raw) == nrow(fit_data),
    "Domestic cultivation counts differ from the fitted data" =
      all(group_raw$planting.China.tdwg3 == fit_data$cul.China),
    "Foreign cultivation counts differ from the fitted data" =
      all(group_raw$WorldCuP.n.tdwg3 == fit_data$cul.foreign)
  )

  message("g-computation: ", group)
  set.seed(base_seed + i)

  # ---- 3.2 Posterior draws ----
  # fixef(summary = FALSE) returns a draws x coefficients matrix with the
  # chains stacked; as_draws_matrix() uses the same draw order, so draw_index
  # selects matching rows from both. Evenly spaced indices across the stacked
  # chains take roughly equal numbers of draws from every chain. unique()
  # removes indices duplicated by rounding, so n_draws can be slightly smaller
  # than max_draws.
  coef_draws <- fixef(fit, summary = FALSE)  # rows = draws, columns = fixed effects
  draw_index <- unique(round(seq(
    1, nrow(coef_draws),
    length.out = min(max_draws, nrow(coef_draws))
  )))
  coef_draws <- coef_draws[draw_index, , drop = FALSE]

  # Auxiliary (non-regression) parameters
  #   phi_culChina, phi_culforeign : beta-binomial precision; larger values
  #                                  mean less overdispersion
  #   shape_natextent              : NB shape; larger values mean less
  #                                  overdispersion (Poisson as shape -> Inf)
  aux_draws <- as_draws_matrix(
    fit,
    variable = c("phi_culChina", "phi_culforeign", "shape_natextent")
  )
  phi_china   <- as.numeric(aux_draws[draw_index, "phi_culChina"])
  phi_foreign <- as.numeric(aux_draws[draw_index, "phi_culforeign"])
  nb_shape    <- as.numeric(aux_draws[draw_index, "shape_natextent"])
  n_draws     <- nrow(coef_draws)

  # ---- 3.3 Recover the standardization used in fitting ----
  # The model was fitted on z = (log1p(v) - center) / scale, computed within
  # the group. brms does not store center and scale, so they are recomputed
  # from the same values: native range from the raw data, cultivation counts
  # from fit$data. sd() uses the n - 1 denominator, as scale() does.
  # Simulated counts are later mapped with these same constants, so that they
  # are on the scale of the fitted coefficients.
  native_center  <- mean(log1p(native_range_raw))
  native_scale   <- sd(log1p(native_range_raw))
  china_center   <- mean(log1p(fit_data$cul.China))
  china_scale    <- sd(log1p(fit_data$cul.China))
  foreign_center <- mean(log1p(fit_data$cul.foreign))
  foreign_scale  <- sd(log1p(fit_data$cul.foreign))

  # Stop if the recomputed z-scores differ from the fitted predictors, e.g.
  # because the predictors were standardized on another subset or scale
  stopifnot(
    "Native range standardization differs from model fitting" = max(abs(
      (log1p(native_range_raw) - native_center) / native_scale -
        fit_data$native.range.scaled
    )) < 1e-10,
    "Domestic cultivation standardization differs from model fitting" = max(abs(
      (log1p(fit_data$cul.China) - china_center) / china_scale -
        fit_data$cul.China.scaled
    )) < 1e-10,
    "Foreign cultivation standardization differs from model fitting" = max(abs(
      (log1p(fit_data$cul.foreign) - foreign_center) / foreign_scale -
        fit_data$cul.foreign.scaled
    )) < 1e-10
  )

  # ---- 3.4 Contrast levels and binomial trial numbers ----
  # Q25 and Q75 are taken on the raw count scale (quantile() type 7, i.e.
  # linear interpolation) and only then log1p-transformed and standardized.
  native_q25_q75 <- quantile(native_range_raw, c(0.25, 0.75), names = FALSE)

  # Number of beta-binomial trials, i.e. the number of regions available for
  # cultivation in China and abroad (fixed in fitting, e.g. 10 and 359).
  # A single value is used in the simulations, so it must be constant across
  # species.
  trials_china   <- unique(fit_data$trials.China)
  trials_foreign <- unique(fit_data$trials.foreign)
  stopifnot(
    "trials.China is not constant"   = length(trials_china) == 1L,
    "trials.foreign is not constant" = length(trials_foreign) == 1L
  )

  # ---- 3.5 Common random numbers ----
  # One set of uniforms for the Beta step and one for the Binomial step of
  # each mediator (n_draws * n_sims values each). They are turned into draws
  # by inverse-CDF sampling (qbeta, qbinom) in Sections 4 and 5.
  #   - The China and foreign sets are drawn independently. This matches the
  #     model, in which the two cultivation equations have independent errors
  #     given their predictors.
  #   - Within a mediator, every scenario re-uses the same set (the coupling
  #     described in the header); the null check below re-uses them as well.
  # The order of these four lines fixes the random-number stream; reordering
  # them changes the results.
  u_china_beta    <- runif(n_draws * n_sims)
  u_china_binom   <- runif(n_draws * n_sims)
  u_foreign_beta  <- runif(n_draws * n_sims)
  u_foreign_binom <- runif(n_draws * n_sims)

  # Additivity errors of the actual contrast, stored so that the null check
  # can report the larger of the two errors (Section 7.1)
  additivity_error_main <- setNames(rep(NA_real_, length(outcomes)), outcomes)

  # ---- Contrast loop ----
  # Q25_to_Q75      : the actual contrast (x0 = Q25, x1 = Q75)
  # No_change_check : a null check with x1 = x0. Every scenario then receives
  #                   identical inputs and random numbers, so every path
  #                   effect must be zero (exactly, since the computations
  #                   are identical). A non-zero value would reveal an error
  #                   in how the scenarios are wired.
  # Both runs go through identical code. Q25_to_Q75 must run first, because
  # the check reads its additivity errors.
  for (comparison in c("Q25_to_Q75", "No_change_check")) {
    native_low  <- native_q25_q75[1]
    native_high <- if (comparison == "No_change_check") {
      native_low
    } else {
      native_q25_q75[2]
    }

    # Standardized X at the two levels; the scenario tables refer to these
    # values by the names "low" and "high"
    x_scaled <- c(
      low  = (log1p(native_low)  - native_center) / native_scale,
      high = (log1p(native_high) - native_center) / native_scale
    )

    # ========================================================================
    # 4. Simulate domestic cultivation C at x0 and x1
    #    china_sims$low = C(x0), china_sims$high = C(x1)
    #    logit(p_C) = c0 + c1 * X
    #    latent ~ Beta(p_C * phi, (1 - p_C) * phi)   (mean p_C, precision phi)
    #    count  ~ Binomial(trials, latent)
    # ========================================================================
    # Vector recycling: p and phi_china have length n_draws, whereas the
    # uniforms have length n_draws * n_sims. R recycles p and phi along the
    # uniforms, so element k uses posterior draw ((k - 1) %% n_draws) + 1.
    # matrix() fills by column with nrow = n_draws, which puts element k in
    # that same row. Hence row r always corresponds to posterior draw r, and
    # each column is one mediator simulation.
    #
    # Inverse-CDF sampling: for a fixed uniform, qbeta() increases with p_C
    # and qbinom() increases with the latent probability. Whichever level of
    # X gives the higher p_C therefore never gives a smaller count for the
    # same uniform, which is what makes C(x0) and C(x1) rank-correlated.
    china_sims <- list()
    for (level in c("low", "high")) {
      # Mean cultivation probability for each posterior draw (length n_draws)
      p <- plogis(
        coef_draws[, "culChina_Intercept"] +
          coef_draws[, "culChina_native.range.scaled"] * x_scaled[[level]]
      )
      # Keep both Beta shape parameters strictly positive
      p <- pmin(pmax(p, 1e-9), 1 - 1e-9)
      latent_p <- qbeta(u_china_beta, p * phi_china, (1 - p) * phi_china)
      counts <- matrix(
        qbinom(u_china_binom, size = trials_china, prob = latent_p),
        nrow = n_draws, ncol = n_sims    # rows = draws, columns = simulations
      )
      # Map simulated counts (0 to trials_china) onto the model scale
      china_sims[[level]] <- (log1p(counts) - china_center) / china_scale
    }

    # ========================================================================
    # 5. Simulate foreign cultivation F for the three F scenarios
    #    logit(p_F) = f0 + f1 * X + f2 * C   (scenarios defined in Section 2)
    #    Sampling as in Section 4, using the foreign uniforms
    # ========================================================================
    # The linear predictor is built as an n_draws x n_sims matrix:
    #   - matrix(f0 + f1 * x, ...) copies the per-draw vector into every
    #     column, because X is the same in all simulations;
    #   - sweep(china_sim, 1, f2, "*") multiplies row r of the simulated C by
    #     the coefficient of posterior draw r, because C varies by simulation.
    # p is now a full matrix, so only phi_foreign (length n_draws) is
    # recycled, again along the rows.
    foreign_sims <- list()
    for (j in seq_len(nrow(foreign_scenarios))) {
      scenario  <- foreign_scenarios[j, ]
      x         <- x_scaled[[scenario$x]]          # standardized X in F
      china_sim <- china_sims[[scenario$china]]    # simulated C entering F

      eta <- matrix(
        coef_draws[, "culforeign_Intercept"] +
          coef_draws[, "culforeign_native.range.scaled"] * x,
        nrow = n_draws, ncol = n_sims
      ) +
        sweep(china_sim, 1, coef_draws[, "culforeign_cul.China.scaled"], "*")
      p <- pmin(pmax(plogis(eta), 1e-9), 1 - 1e-9)
      latent_p <- qbeta(u_foreign_beta, p * phi_foreign, (1 - p) * phi_foreign)
      counts <- matrix(
        qbinom(u_foreign_binom, size = trials_foreign, prob = latent_p),
        nrow = n_draws, ncol = n_sims
      )
      # Map simulated counts (0 to trials_foreign) onto the model scale
      foreign_sims[[scenario$name]] <- (log1p(counts) - foreign_center) / foreign_scale
    }

    # ========================================================================
    # 6. Outcome predictions for the five Y scenarios
    #    The three estimands share the same mediator simulations.
    #    Linear predictors of the two hurdle parts:
    #      log(mu)   = a0 + a1 * X + a2 * C + a3 * F
    #      logit(hu) = h0 + h1 * X + h2 * C + h3 * F
    # ========================================================================
    expected_y <- list(incidence = list(), extent = list(), success = list())
    for (j in seq_len(nrow(outcome_scenarios))) {
      scenario    <- outcome_scenarios[j, ]
      x           <- x_scaled[[scenario$x]]              # X in Y
      china_sim   <- china_sims[[scenario$china]]        # C in Y
      foreign_sim <- foreign_sims[[scenario$foreign]]    # F in Y

      # Both linear predictors are built as in Section 5: the X term is
      # copied across columns, and the C and F terms are scaled row-wise by
      # the coefficients of each posterior draw.
      eta <- list()
      for (part in names(hurdle_prefix)) {
        prefix <- hurdle_prefix[[part]]
        eta[[part]] <- matrix(
          coef_draws[, paste0(prefix, "Intercept")] +
            coef_draws[, paste0(prefix, "native.range.scaled")] * x,
          nrow = n_draws, ncol = n_sims
        ) +
          sweep(china_sim, 1,
                coef_draws[, paste0(prefix, "cul.China.scaled")], "*") +
          sweep(foreign_sim, 1,
                coef_draws[, paste0(prefix, "cul.foreign.scaled")], "*")
      }

      mu <- exp(eta$mu)     # mean of the untruncated NB distribution
      hu <- plogis(eta$hu)  # probability of a zero

      # Probability of a positive count under the untruncated NB:
      #   1 - P_NB(0) = 1 - (1 + mu / shape)^(-shape)
      #               = -expm1(-shape * log1p(mu / shape))
      # The expm1/log1p form avoids loss of precision when mu is small
      # (P_NB(0) close to 1). nb_shape (length n_draws) is recycled along
      # the rows, as in Section 4.
      p_positive <- -expm1(-nb_shape * log1p(mu / nb_shape))

      # Estimands for every draw x simulation cell
      incidence <- 1 - hu                       # Pr(Y > 0)
      extent    <- mu / p_positive              # E[Y | Y > 0]
      success   <- incidence * mu / p_positive  # E[Y]; same order of operations
                                                # as (1 - hu) * mu / (1 - P0)

      # Stop on overflow (very large mu) or 0/0 (mu underflowing to zero)
      stopifnot(
        "Non-finite predicted values" =
          all(is.finite(incidence), is.finite(extent), is.finite(success))
      )

      # Average over mediator simulations (columns). This integrates over the
      # simulated distribution of C and F and gives one expected value per
      # posterior draw (a vector of length n_draws).
      expected_y$incidence[[scenario$name]] <- rowMeans(incidence)
      expected_y$extent[[scenario$name]]    <- rowMeans(extent)
      expected_y$success[[scenario$name]]   <- rowMeans(success)
    }

    # ========================================================================
    # 7. Path decomposition
    #    Each outcome is decomposed and normalized on its own; the outcomes
    #    are not expressed relative to one another.
    # ========================================================================
    for (outcome in outcomes) {
      y <- expected_y[[outcome]]  # named list: one vector per scenario

      # Path effects for every posterior draw (n_draws x 4 matrix):
      # differences between adjacent scenarios
      path_draws <- cbind(
        Direct      = y$Direct  - y$Baseline,
        Via_China   = y$China   - y$Direct,
        Via_Foreign = y$Foreign - y$China,
        Via_Both    = y$Both    - y$Foreign
      )
      total_draws <- y$Both - y$Baseline

      # The four differences telescope to the total effect, so they must add
      # up to it exactly, apart from floating-point rounding
      additivity_error <- max(abs(rowSums(path_draws) - total_draws))
      stopifnot(
        "The four paths do not sum to the total effect" = additivity_error < 1e-8
      )
      key <- paste(group, outcome, sep = "_")

      # ---- 7.1 Null check: with X unchanged, every path must be zero ----
      # Only the check summary is stored; `next` then skips Section 7.2.
      if (comparison == "No_change_check") {
        null_error <- max(abs(path_draws))
        stopifnot(
          "Path effects are non-zero when X is unchanged" = null_error < 1e-10
        )
        check_tables[[key]] <- tibble(
          Group               = group,
          Outcome             = outcome,
          Null_effect_max_abs = null_error,
          Additivity_max_abs  = max(additivity_error_main[[outcome]],
                                    additivity_error)
        )
        next
      }

      # ---- 7.2 Actual contrast: posterior means and relative contributions ----
      # Effects are posterior means of the per-draw differences; no credible
      # intervals are computed. Shares are computed from these means (a ratio
      # of means) rather than averaged over draws, which would be unstable in
      # draws where the total effect is close to zero.
      additivity_error_main[[outcome]] <- additivity_error
      path_effects <- colMeans(path_draws)  # named vector of the four paths
      total_effect <- mean(total_draws)

      # Absolute share = |effect| / sum(|effects|). Always between 0% and
      # 100%, and the four paths sum to 100%. Shows how large each path is
      # relative to all paths combined, regardless of sign.
      # NA if every path is exactly zero.
      share_abs <- if (sum(abs(path_effects)) > 0) {
        100 * abs(path_effects) / sum(abs(path_effects))
      } else {
        rep(NA_real_, 4)
      }
      # Signed share = effect / total effect. The four paths also sum to
      # 100%, but a single path can be negative or exceed 100% when paths act
      # in opposite directions. NA if the total effect is (near) zero.
      share_signed <- if (abs(total_effect) > 1e-10) {
        100 * path_effects / total_effect
      } else {
        rep(NA_real_, 4)
      }

      # One row per path. Native_low and Native_high are on the raw scale;
      # N_MC is the number of mediator simulations per draw (n_sims).
      effect_tables[[key]] <- tibble(
        Group               = group,
        Outcome             = outcome,
        Path                = names(path_effects),
        Effect              = unname(path_effects),
        Relative_abs_pct    = unname(share_abs),
        Relative_signed_pct = unname(share_signed),
        Total_effect        = total_effect,
        Unit                = outcome_units[[outcome]],
        Native_low          = native_low,
        Native_high         = native_high,
        N_species           = nrow(fit_data),
        N_draws             = n_draws,
        N_MC                = n_sims
      )

      # Posterior mean of each of the five scenarios, useful for checking the
      # levels behind each difference
      scenario_tables[[key]] <- tibble(
        Group          = group,
        Outcome        = outcome,
        Scenario       = names(y),
        Expected_value = unname(map_dbl(y, mean))
      )
    }
  }
}

# ============================================================================
# 8. Combine, check and export
# ============================================================================

# Stack the per-group, per-outcome tables into three long tables
gcomp_effects   <- bind_rows(effect_tables)
gcomp_scenarios <- bind_rows(scenario_tables)
gcomp_checks    <- bind_rows(check_tables)

# Expected sizes: 4 groups x 3 outcomes x 4 paths = 48 effect rows;
# 4 groups x 3 outcomes = 12 check rows
stopifnot(
  "Effect table does not have 48 rows" = nrow(gcomp_effects) == 48L,
  "Check table does not have 12 rows"  = nrow(gcomp_checks) == 12L
)

# Absolute shares must sum to 100% within each group x outcome. Rows with NA
# shares (all paths zero) are dropped first, mirroring the default NA removal
# of aggregate() in the earlier version of this script.
share_totals <- gcomp_effects %>%
  filter(!is.na(Relative_abs_pct)) %>%
  group_by(Group, Outcome) %>%
  summarise(share_sum = sum(Relative_abs_pct), .groups = "drop")
stopifnot(
  "Absolute shares do not sum to 100%" =
    all(abs(share_totals$share_sum - 100) < 1e-8)
)

# write_csv() writes the same values as write.csv(), but does not quote
# strings and keeps enough digits to read the numbers back exactly. Use
# write.csv(..., row.names = FALSE) if byte-identical files are needed.
write_csv(gcomp_effects,
          file.path(result_dir, "gcomp_three_outcomes_effects.csv"))
# Optional outputs: scenario means and check summaries
# write_csv(gcomp_scenarios,
#           file.path(result_dir, "gcomp_three_outcomes_scenarios.csv"))
# write_csv(gcomp_checks,
#           file.path(result_dir, "gcomp_three_outcomes_checks.csv"))

message("Results directory: ", normalizePath(result_dir, winslash = "/"))

# ============================================================================
# 9. R^2 for each species group: the two cultivation equations and each
#    hurdle part
#    Requires the following externally defined functions (run the script
#    that defines them before this section):
#      residual_R2()   : takes observed values and a draws x species matrix
#                        of expected values (as called below)
#      summarise_R2()  : summarizes the output of residual_R2() for one
#                        model component
#      get_hurdle_R2() : R^2 of the hurdle parts of the naturalization model
# ============================================================================
# pred: posterior draws × observations
residual_R2 <- function(y, pred) {
  pred <- as.matrix(pred)
  stopifnot(length(y) == ncol(pred), length(y) >= 2L,
            all(is.finite(y)), all(is.finite(pred)))
  
  var.fit <- apply(pred, 1, var)
  var.res <- apply(sweep(pred, 2, y, FUN = "-"), 1, var)
  denominator <- var.fit + var.res
  
  out <- rep(NA_real_, length(denominator))
  valid <- denominator > 0
  out[valid] <- var.fit[valid] / denominator[valid]
  out
}

summarise_R2 <- function(x, component, n) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(tibble(component = component, n = n,
                  R2_median = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_))
  }
  tibble(component = component, n = n,
         R2_median = median(x),
         CI_lower = unname(quantile(x, 0.025)),
         CI_upper = unname(quantile(x, 0.975)))
}

get_hurdle_R2 <- function(fit) {
  required <- c(
    "nat.extent", "cul.China", "cul.foreign",
    "native.range.scaled", "cul.China.scaled",
    "cul.foreign.scaled", "trials.China", "trials.foreign"
  )
  
  dat <- fit$data
  dat <- dat[complete.cases(dat[, required, drop = FALSE]), , drop = FALSE]
  
  y <- dat$nat.extent
  incidence <- as.integer(y > 0)
  positive <- y > 0
  
  hu    <- posterior_epred(fit, newdata = dat, resp = "natextent", dpar = "hu")
  mu    <- posterior_epred(fit, newdata = dat, resp = "natextent", dpar = "mu")
  shape <- posterior_epred(fit, newdata = dat, resp = "natextent", dpar = "shape")
  
  p.incidence <- 1 - hu
  prob.positive.NB <- -expm1(-shape * log1p(mu / shape))
  mean.extent <- mu / prob.positive.NB
  mean.success <- p.incidence * mean.extent
  
  bind_rows(
    summarise_R2(residual_R2(incidence, p.incidence),
                 "Hurdle_incidence", length(incidence)),
    summarise_R2(residual_R2(y[positive], mean.extent[, positive, drop = FALSE]),
                 "Hurdle_extent", sum(positive)),
    summarise_R2(residual_R2(y, mean.success),
                 "Hurdle_success", length(y))
  )
}
# R^2 of the two cultivation equations.
# posterior_epred() returns, for each draw and species, the expected count on
# the response scale (trials x p), not the logit-scale linear predictor, so
# the predictions can be compared directly with the observed counts. It uses
# all posterior draws, not the thinned set from Section 3.2.
get_cultivation_R2 <- function(fit) {
  # All variables of the joint model, including the trial columns that
  # newdata requires. brms already drops incomplete rows when fitting, so this
  # is a safeguard that keeps both equations on the same species.
  required_cols <- c(
    "nat.extent", "cul.China", "cul.foreign",
    "native.range.scaled", "cul.China.scaled", "cul.foreign.scaled",
    "trials.China", "trials.foreign"
  )
  r2_data <- fit$data %>%
    drop_na(all_of(required_cols))

  # draws x species matrices of expected cultivation counts
  pred_china   <- posterior_epred(fit, newdata = r2_data, resp = "culChina")
  pred_foreign <- posterior_epred(fit, newdata = r2_data, resp = "culforeign")

  # One summary per cultivation equation, stacked into a single table
  bind_rows(
    summarise_R2(
      residual_R2(r2_data$cul.China, pred_china),
      component = "Cultivation_China", n = nrow(r2_data)
    ),
    summarise_R2(
      residual_R2(r2_data$cul.foreign, pred_foreign),
      component = "Cultivation_foreign", n = nrow(r2_data)
    )
  )
}

# Apply to every group: imap() passes each model together with its list name,
# which is added as the first column (Group)
r2_results <- models %>%
  imap(function(fit, group) {
    bind_rows(get_cultivation_R2(fit), get_hurdle_R2(fit)) %>%
      mutate(Group = group, .before = 1)
  }) %>%
  list_rbind()

print(r2_results, n = Inf)  # print all rows
write_csv(r2_results,
          file.path(result_dir, "R2_all_equations_and_hurdle_parts.csv"))