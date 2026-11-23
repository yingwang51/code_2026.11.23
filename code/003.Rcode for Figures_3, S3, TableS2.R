# ============================================================================
# Two-part Bayesian SEM for naturalization
# No user-defined functions: the analysis is written step by step.
# ============================================================================

cat("\014")
rm(list = ls())

library(brms)
library(dplyr)
library(posterior)
library(readr)
library(tidyr)
library(ggplot2)

setwd("E:/中国本地植物归化/sem分析")

data.file <- "./20250530.native_plant.matched.cultivated_plant_DO11.csv"
result.dir <- "./2026.1104.results"
dir.create(result.dir, recursive = TRUE, showWarnings = FALSE)

# Number of possible TDWG level-3 cultivation regions.
n.China <- 10L
n.foreign.total <- 359L
seed.use <- 2024

# ============================================================================
# 1. Data preparation
# ============================================================================
# Remove data without life forms
native.species <- read_csv(data.file, show_col_types = FALSE) %>%
  filter(!is.na(life.form.integrated)) 

native.species00 <- native.species %>%
  transmute(
    life.form = life.form.integrated,
    native.range = as.integer(native.global.tdwg3),
    cul.China = as.integer(planting.China.tdwg3),
    cul.foreign = as.integer(WorldCuP.n.tdwg3),
    nat.incidence = as.integer(nat.incidence),
    nat.extent = as.integer(nat.extent),
    trials.China = n.China,
    trials.foreign = n.foreign.total
  )

# ============================================================================
# 2. Step 1: naturalization incidence
#    All species; Bernoulli response
# ============================================================================

# All species
native.incidence00 <- native.species00 %>%
  mutate(
    native.range.scaled = as.numeric(scale(log1p(native.range))),
    cul.China.scaled = as.numeric(scale(log1p(cul.China))),
    cul.foreign.scaled = as.numeric(scale(log1p(cul.foreign)))
  )

# Annual herbs
native.incidence01 <- native.species00 %>%
  filter(life.form == "annual herb") %>%
  mutate(
    native.range.scaled = as.numeric(scale(log1p(native.range))),
    cul.China.scaled = as.numeric(scale(log1p(cul.China))),
    cul.foreign.scaled = as.numeric(scale(log1p(cul.foreign)))
  )

# Perennial herbs
native.incidence02 <- native.species00 %>%
  filter(life.form == "perennial herb") %>%
  mutate(
    native.range.scaled = as.numeric(scale(log1p(native.range))),
    cul.China.scaled = as.numeric(scale(log1p(cul.China))),
    cul.foreign.scaled = as.numeric(scale(log1p(cul.foreign)))
  )

# Woody plants
native.incidence03 <- native.species00 %>%
  filter(life.form == "woody") %>%
  mutate(
    native.range.scaled = as.numeric(scale(log1p(native.range))),
    cul.China.scaled = as.numeric(scale(log1p(cul.China))),
    cul.foreign.scaled = as.numeric(scale(log1p(cul.foreign)))
  )

# Full-species reference data reused by both stages in g-computation.
reference.list <- list(
  All = native.incidence00,
  Annual = native.incidence01,
  Perennial = native.incidence02,
  Woody = native.incidence03
)


# ============================================================================
# naturalization success BSEM formulas
formula.success.cul.China <- bf(
  cul.China | trials(trials.China) ~ native.range.scaled,
  family = beta_binomial(link = "logit")
)

formula.success.cul.foreign <- bf(
  cul.foreign | trials(trials.foreign) ~
    cul.China.scaled + native.range.scaled,
  family = beta_binomial(link = "logit")
)


# Naturalization: hurdle negative binomial model
formula.nat.success <- bf(
  # 正计数部分：零截断负二项分布
  nat.extent ~
    cul.foreign.scaled + cul.China.scaled + native.range.scaled,
  
  # 零值部分：预测 P(nat.extent = 0)
  hu ~
    cul.foreign.scaled + cul.China.scaled + native.range.scaled,
  
  family = hurdle_negbinomial(
    link = "log",
    link_hu = "logit"
  )
)
total.cores <- parallel::detectCores(logical = TRUE)

# Model 00: all naturalized species
model.success00 <- brm(
  formula.success.cul.China +
    formula.success.cul.foreign +
    formula.nat.success +
    set_rescor(FALSE),
  data = native.incidence00,
  chains = 4,
  cores = total.cores/4,
  iter = 2000,
  seed = seed.use
)

# Model 01: naturalized annual herbs
model.success01 <- brm(
  formula.success.cul.China +
    formula.success.cul.foreign +
    formula.nat.success +
    set_rescor(FALSE),
  data = native.incidence01,
  chains = 4,
  cores = total.cores/4,
  iter = 2000,
  seed = seed.use
)

# Model 02: naturalized perennial herbs
model.success02 <- brm(
  formula.success.cul.China +
    formula.success.cul.foreign +
    formula.nat.success +
    set_rescor(FALSE),
  data = native.incidence02,
  chains = 4,
  cores = total.cores/4,
  iter = 2000,
  seed = seed.use
)

# Model 03: naturalized woody plants
model.success03 <- brm(
  formula.success.cul.China +
    formula.success.cul.foreign +
    formula.nat.success +
    set_rescor(FALSE),
  data = native.incidence03,
  chains = 4,
  cores = total.cores/4,
  iter = 2000,
  seed = seed.use
)

summary(model.success00)
summary(model.success01)
summary(model.success02)
summary(model.success03)

save(
  model.success00,
  model.success01,
  model.success02,
  model.success03,
  file = file.path(result.dir, "models_naturalization_success.RData")
)
