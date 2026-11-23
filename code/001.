# ==============================================================================
# Manuscript:   Cultivation facilitates global naturalization
# Description:  Collinearity diagnostics for the binomial and zero-truncated
#               negative binomial components of the naturalization model.
# Author:       dbc and wy
# Output:       results/Table S3.csv
# ==============================================================================

# 1. Environment Setup ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr,   # Read and write CSV files
  dplyr,   # Data manipulation
  tibble,  # Convert row names to a column
  VGAM,    # Zero-truncated negative binomial models
  car      # Generalized variance inflation factors
)

rm(list = ls())
cat("\014")

# --- Project Paths ---
# Change project_dir if the data are stored elsewhere.
project_dir <- "your_work_path"
data_file <- file.path(
  project_dir, "20250530.native_plant.matched.cultivated_plant_DO11.csv"
)
result_dir <- file.path(project_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Data Loading & Preprocessing ----------------------------------------------
native_raw <- read_csv(data_file, show_col_types = FALSE)

# --- Select Variables and Set the Reference Life Form ---
native_clean <- native_raw %>%
  select(
    taxon_name,
    nat.extent,
    planting.China.tdwg3,
    WorldCuP.n.tdwg3.NEW,
    native.global.tdwg3,
    life.form.integrated
  ) %>%
  filter(!is.na(life.form.integrated)) %>%
  mutate(
    life.form.integrated = factor(
      life.form.integrated,
      levels = c("woody", "annual herb", "perennial herb")
    )
  )

# --- Transform and Standardize Predictors ---
# Rationale: Log-transformation reduces skewness; z-scaling puts predictors
#            on a common scale. Woody plants are the reference life form.
native_scaled <- native_clean %>%
  mutate(
    scale.WorldCuP.n.tdwg3 = as.numeric(scale(log(WorldCuP.n.tdwg3.NEW + 1))),
    scale.native.global.tdwg3 = as.numeric(scale(log(native.global.tdwg3 + 1))),
    scale.planting.China.tdwg3 = as.numeric(scale(log(planting.China.tdwg3 + 1)))
  )

# Check the processed data and raw predictor ranges.
str(native_scaled)
range(native_scaled$WorldCuP.n.tdwg3.NEW)
range(native_scaled$native.global.tdwg3)
range(native_scaled$planting.China.tdwg3)

# 3. Statistical Modeling ------------------------------------------------------

# --- Model Formula: Main Effects and Life-Form Interactions ---
model_formula <- nat.extent ~
  scale.native.global.tdwg3 +
  scale.WorldCuP.n.tdwg3 +
  scale.planting.China.tdwg3 +
  life.form.integrated +
  scale.native.global.tdwg3:life.form.integrated +
  scale.WorldCuP.n.tdwg3:life.form.integrated +
  scale.planting.China.tdwg3:life.form.integrated

# --- Binomial Model: Naturalization Incidence ---
zero_glm <- glm(
  model_formula,
  data = native_scaled %>% mutate(nat.extent = as.numeric(nat.extent > 0)),
  family = binomial
)
summary(zero_glm)

# --- Zero-Truncated Negative Binomial Model: Naturalization Extent ---
data_nonzero <- subset(native_scaled, nat.extent > 0)
count_vglm <- vglm(
  model_formula,
  data = data_nonzero,
  family = posnegbinomial
)
summary(count_vglm)

# 4. Collinearity Diagnostics --------------------------------------------------

# --- Binomial Component: Generalized VIF ---
vif_zero <- car::vif(zero_glm, type = "predictor")

df_vif_zero <- as.data.frame(vif_zero) %>%
  rownames_to_column(var = "Term") %>%
  rename(GVIF_adj = `GVIF^(1/(2*Df))`) %>%
  mutate(Part = "Zero")

# --- Function: Calculate Column-Wise VIF for the Count Component ---
# Rationale: Regress each retained design-matrix column on the remaining columns,
#            then calculate VIF = 1 / (1 - R2).
calc_vif_vglm <- function(model) {
  X <- model.matrix(model)
  # Remove columns containing only zeros.
  X <- X[, !apply(X, 2, function(x) all(x == 0)), drop = FALSE]

  vif_values <- numeric(ncol(X))
  names(vif_values) <- colnames(X)

  for (i in seq_len(ncol(X))) {
    y <- X[, i]
    x <- X[, -i, drop = FALSE]
    r2 <- summary(lm(y ~ x))$r.squared
    vif_values[i] <- 1 / (1 - r2)
  }

  return(vif_values)
}

# --- Count Component: Calculate and Organize VIF ---
vif_count <- calc_vif_vglm(count_vglm)

df_vif_count <- data.frame(
  Term = names(vif_count),
  GVIF = vif_count,
  Df = 1,
  GVIF_adj = sqrt(vif_count),
  Part = "Count"
)

# 5. Combine and Export Results ------------------------------------------------
df_vif <- bind_rows(df_vif_count, df_vif_zero)
print(df_vif)

write_csv(df_vif, file.path(result_dir, "Table_S3.csv"))

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
