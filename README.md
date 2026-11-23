## README: Data and R Code for "Associations of native range size and cultivation with global naturalization vary among growth habits of plants native to China"

---

### Data

**1. 20250530.native_plant.matched.cultivated_plant_DO11.csv**

**Description**: This is the core dataset for analyzing the global naturalization of plants native to China. It integrates taxonomic identifiers, growth habits, naturalization status and extent, native range size, and cultivation status and extent within and outside China.

*Taxonomic identifiers*

- `accepted_plant_name_id`: The identifier of the accepted name in WCVP
- `taxon_name`: The accepted name (Latin name) of the species
- `taxon_authors`: The nomenclatural authorship of the accepted scientific name

*Growth habits, cultivation status, and naturalization status*

- `life.form.integrated`: Integrated growth habit of the species (annual herb, perennial herb, or woody)
- `cultivation.lin2019`: Cultivation status determined according to *Catalogue of Cultivated Plants in China* (yes or no)
- `nat.extent`: Naturalization extent (number of regions where plants are naturalized)
- `nat.incidence`: Binary naturalization status (0 = not naturalized; 1 = naturalized)

*Native range size and cultivation extent*

- `native.global.tdwg3`: Native range size (number of TDWG level-3 regions where the species is native)
- `planting.China.tdwg3`: Cultivation within China, measured as the number of TDWG level-3 regions within China where the species is recorded as cultivated
- `WorldCuP.n.tdwg3`: Cultivation outside China, measured as the number of TDWG level-3 regions outside China where the species is recorded as cultivated

**2. NativeChinesePlantsCultivationExtentWorldCuP.csv**

**Description**: This dataset summarizes WorldCuP cultivation records for 25,808 taxa native to China, of which 9,670 have at least one WorldCuP list record. Cultivation extent outside China is recorded as the number of TDWG level-3 regions.

- `taxon_name`: The accepted name (Latin name) of the species
- `taxon_authors`: The nomenclatural authorship of the accepted scientific name
- `accepted_plant_name_id`: The unique identifier for the accepted name
- `number_records`: The number of unique WorldCuP lists in which each taxon is present
- `number_tdwg3`: The total number of unique TDWG level-3 regions outside of China where the taxon is recorded as cultivated according to WorldCuP

**3. References_WorldCuPSources.csv**

**Description**: This file contains the comprehensive list of 365 bibliographic references and data sources used to compile the cultivation extent information for Chinese native plants within the World Checklist of Cultivated Plants (WorldCuP) database.

---

### Code

#### 1. Rcode for Table_S3.R

**Purpose**: Collinearity diagnostics for the binomial and zero-truncated negative binomial components of the naturalization model.

**Key analyses**:

- Log-transform and z-standardize native range size and cultivation extent within and outside China, with woody plants as the reference growth habit
- Fit binomial and zero-truncated negative binomial models including growth habit and its interactions with the three continuous predictors
- Calculate generalized variance inflation factors (GVIFs) for the binomial model using `car::vif()` and column-wise VIFs for the count model using auxiliary regressions on its design matrix
- Combine the diagnostics, degrees of freedom, and adjusted values in a single table, identifying the model component for each entry

**Output**: Table S3 (collinearity diagnostics for both model components; `Table_S3.csv`)

---

#### 2. Rcode for Figures_3, S3, TableS2.R

**Purpose**: Bayesian structural equation modeling of associations among native range size, cultivation extent within and outside China, and naturalization across all species and the three growth habit groups.

**Key analyses**:

- Fit a multivariate Bayesian SEM separately for all species, annual herbs, perennial herbs, and woody plants, using beta-binomial cultivation responses and a hurdle negative binomial naturalization response
- Derive three views from each fitted model: naturalization incidence from the sign-reversed hurdle coefficients, extent from the positive-count component, and overall success by combining the hurdle and count components
- Estimate overall success effects as average direct slopes of log expected naturalization extent across observed predictor profiles
- Visualize path estimates on their corresponding effect scales, with posterior medians, 95% credible intervals, and directional probabilities

**Outputs**:

- Figures 3 and S3: BSEM path diagrams for naturalization incidence, extent, and overall success, each covering all species and the three growth habit groups
  
- Table S2: Posterior summaries for 72 plotted paths (3 outcomes × 4 groups × 6 paths), including effect estimates, 95% credible intervals, directional probabilities, and effect scales
  

---

#### 3. Rcode for Figures_2, S1 and Tables_1, S1.R

**Purpose**: Hurdle model analysis examining how native range size, cultivation extent within and outside China, and growth habit predict naturalization success, incidence, and extent.

**Key analyses**:

- Fit a negative binomial hurdle model with growth habit interactions in both the binomial and zero-truncated negative binomial components
- Fit component models separately for detailed inspection:
  - Binomial model: probability of naturalization (yes/no)
  - Zero-truncated negative binomial: naturalization extent (number of TDWG3 regions, conditional on occurrence)
- Generate model predictions across focal variable gradients
- Create visualization panels with inset zoom plots for low-range details
- Predict growth habit effects while holding continuous predictors at their mean values, using positive-count data for the count component
- Perform 14 likelihood ratio tests of main effects and growth habit interactions across the two model components; test main effects after removing interactions from the component being examined

**Outputs**:

- Figure 2: Predicted relationships of native range size and cultivation extent with naturalization success, incidence, and extent, including the combined nine-panel visualization
- Figure S1: Predicted growth habit effects on naturalization success, incidence, and extent at mean continuous predictor values
- Table 1: Likelihood ratio tests for main effects and growth habit interactions in both hurdle model components
- Table S1: Coefficient estimates, standard errors, test statistics, and p-values from the full hurdle and separate component models

---

#### 4. Rcode for Figure_S2.R

**Purpose**: Visualization of linearized relationships from two-part hurdle models examining predictors of naturalization extent.

**Key analyses**:

- Fit negative binomial hurdle models with growth habit interactions
- Zero component: model probability of naturalization occurrence (binomial regression)
- Count component: model naturalization breadth conditional on occurrence (zero-truncated negative binomial)
- Generate linear predictor plots across focal variables (native range size, cultivation within China, and cultivation outside China)

**Output**: Figure S2 (6-panel plot: 3 predictors × 2 model components)

---

#### 5. Rcode for Figure_S5.R

**Purpose**: Pairwise correlation analysis of native range size, cultivation within China, and cultivation outside China.

**Key analyses**:

- Compute Pearson correlations among log-transformed, z-standardized predictors
- Visualize correlation matrix with scatterplots (lower triangle), correlation coefficients (upper triangle), and marginal distributions (diagonal)
- Test hypothesis that native range size positively correlates with both cultivation within China and cultivation outside China

**Output**: Figure S5 (correlation matrix with regression lines and distribution histograms)
