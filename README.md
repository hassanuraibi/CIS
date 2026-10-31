CIS: Central Information Score for Gamma Regression

R research scripts for Central Information Score (CIS) calibration of the shrinkage parameter in Gamma regression under multicollinearity and contamination. The repository includes a Monte Carlo simulation and an application to the original and jointly perturbed bodyfat data.

Repository files

File

Purpose

cis_core.R

Design standardization, Gamma estimation, screening and refitting, shrinkage solvers, and Excel export utilities.

cis_simulation.R

Defines cis_simulate() for pilot and full simulation studies.

cis_real_data.R

Defines cis_bodyfat() for the fixed bodyfat perturbation study.

These are standalone scripts, not an installable R package. Sourcing the files defines functions; it does not automatically launch an analysis.

Requirements

Install R and the required packages:

install.packages(c("MASS", "openxlsx", "Rmpfr", "TH.data"))

MASS: multivariate normal generation and MVE covariance estimation.

openxlsx: formatted Excel workbooks.

Rmpfr: arbitrary-precision recovery for difficult shrinkage solves. Both analysis entry points check for this package at startup.

TH.data: original bodyfat data; unnecessary for simulation alone or when a compatible data frame is supplied explicitly.

Installing Rmpfr from source can require additional system libraries; use an available binary for your R version where possible. Each analysis saves SessionInfo.txt to record the actual software environment.

Quick start

Place all three scripts in the same folder and set that folder as the R working directory. In RStudio, opening an R project in that folder is convenient. Then run:

source("cis_core.R")
source("cis_simulation.R")
source("cis_real_data.R")

# Small functional run: one setting and two replications.
# This is a software check, not a publishable simulation.
check <- cis_simulate(
  run_mode = "pilot",
  B_override = 2L,
  setting_ids = 1L,
  output_dir = file.path(getwd(), "CIS_Check")
)

# Real-data application: one fixed perturbation, no repeated contamination.
application <- cis_bodyfat(
  output_dir = file.path(getwd(), "bodyfat_Modified")
)

The examples above are usage instructions; they do not imply that the scripts have been validated on every platform. Run the small check before starting the full study.

Simulation study

# Pilot: 20 settings, 50 replications per setting.
pilot <- cis_simulate(
  run_mode = "pilot",
  output_dir = file.path(getwd(), "CIS_Pilot")
)

# Full: 480 settings, 500 replications per setting.
full <- cis_simulate(
  run_mode = "full",
  output_dir = file.path(getwd(), "CIS_Full")
)

The full run attempts 240,000 generated datasets, each evaluated by four methods. High-precision recovery can make some settings slow. Use separate output directories for separate runs: files with the same names can be overwritten. The scripts do not implement automatic resume.

Design implemented in these scripts

Factor

Pilot

Full

Sample size n

50

50, 100

Number of predictors p (excluding intercept)

4, 8

4, 8

Population equicorrelation rho

0.95

0, 0.90, 0.95, 0.99

Gamma dispersion phi

0.50

0.25, 0.50, 1

Coefficient profile

EqualPositive, MixedSigns

EqualPositive, MixedSigns

Conditions

Clean; response contamination at 10% and 20%; joint contamination at 10% and 20%

Same

Replications per setting

50

500

The intercept is 1. EqualPositive uses slopes of 0.5. MixedSigns repeats the slope sequence (0.6, -0.4, 0.3, -0.2) to length p.

Predictors are sampled from a multivariate normal distribution with equicorrelated covariance, then sample-centered and standardized before contamination. Clean responses follow:

$$
Y_i^0\mid x_i^0\sim\operatorname{Gamma}(1/\phi,;\phi\mu_i^0),
\qquad \mu_i^0=\exp{(x_i^0)^\top\beta},
$$

where the Gamma arguments are shape and scale.

For a contamination rate r, ceiling(n * r) distinct rows are sampled. In YOutliers, their responses are multiplied by 8. The code label BadLeverage denotes joint predictor and response contamination: the same selected responses are multiplied by 8, and all standardized predictor columns on those rows are multiplied by 50. It does not denote predictor-only contamination. The intercept remains 1; responses are not regenerated after predictor modification, and predictors are not restandardized.

Replication seeds are:

123L + (SettingID - 1L) * 100000L + replication

To change the number of replications or run a subset:

subset_run <- cis_simulate(
  run_mode = "full",
  B_override = 100L,
  setting_ids = c(1L, 2L, 3L),
  output_dir = file.path(getwd(), "CIS_Subset")
)

Setting identifiers follow the grid order in the script. Settings.rds stores the complete grid before subset filtering; the workbook Settings sheet and the returned settings contain the selected settings. Changing B_override changes replication count, not the factorial grid.

Simulation outputs

Output

Contents

Simulation_Results.xlsx

Settings, Comparison, Performance, Paired, Screening, Failures, NumericalAudit, and Notes worksheets.

Simulation_Figures.pdf

Mean coefficient MSE by population correlation and scenario, where the plotted values are finite and positive.

Settings.rds

Complete ordered settings grid.

Setting_001.rds, etc.

Setting, true coefficients, replication coefficient estimates, losses, screening records, and errors.

SessionInfo.txt

R session and package information.

Simulation coefficient loss is coefficient-averaged MSE, including the intercept:

$$
L_b=\frac{1}{q}\sum_{j=1}^{q}(\widehat\beta_{bj}-\beta_j)^2,
\qquad q=p+1.
$$

Comparison ranks methods using replications successful for all four methods.

Performance reports method-specific successful-replication MSE, squared bias, variance, and failure rate. Empirical variance uses the divisor equal to the number of successful replications.

Paired uses pairwise successful draws. MeanDifference is comparator loss minus CIS loss; positive values favor CIS. Lower and Upper are pointwise Monte Carlo intervals, calculated as the mean difference plus/minus 1.96 Monte Carlo standard errors.

Screening reports retention, tuning, detection, false exclusion, and fallback summaries. Detection and false exclusion use the screening decisions before fallback.

These success subsets can differ. Inspect failure rates alongside MSE; the reported risks are conditional on computational success. The paired intervals are not confidence intervals for regression coefficients.

Bodyfat application

application <- cis_bodyfat(
  contamination_seed = 123L,
  screen_seed = 20261017L,
  output_dir = file.path(getwd(), "bodyfat_Modified")
)

The application loads bodyfat from TH.data: 71 observations, positive response DEXfat, and nine predictors: age, waistcirc, hipcirc, elbowbreadth, kneebreadth, anthro3a, anthro3b, anthro3c, and anthro4.

It selects ceiling(0.10 * 71) = 8 rows once, an actual proportion of approximately 11.27%. On those same rows, all standardized predictors are multiplied by 50 and the response by 8. Original means and standard deviations remain fixed. The modified analysis data must be identified as bodyfat_Modified, not as unaltered observed records.

There are no 500 repeated perturbations and no cross-validation in this application. The 1,000 MVE subsamples are internal covariance-estimation subsets, not contamination replications. The original data are a baseline, not a dataset certified to be free of outliers.

Modified-data fits are evaluated at the original predictor coordinates against the original responses using MAE, RMSE, and mean Gamma unit deviance:

$$
\operatorname{GD}=\frac{2}{n}\sum_{i=1}^{n}
\left[\frac{y_i}{\widehat\mu_i}-1-
\log\left(\frac{y_i}{\widehat\mu_i}\right)\right].
$$

These are response-restoration criteria for one fixed perturbation, not independent test-set performance. No true regression coefficients are assumed. The exported coefficients are in the fixed standardized analysis coordinates; cis_original_beta() is available for conversion to original predictor units when appropriate.

Application outputs

bodyfat_Modified_Results.xlsx: comparison criteria, original/modified coefficient estimates, CIS–RHKB differences, spectral factors, VIFs, screening diagnostics, tuning, status, selected rows, original data, modified analysis coordinates, and notes.

bodyfat_Modified_Figures.pdf: collinearity, robust predictor distances, response-restoration criteria, coefficient proximity, and shrinkage factors.

Fitted_Objects.rds: fits and original/modified analysis matrices and responses.

SessionInfo.txt: software environment.

To supply the original compatible data frame directly:

# original_bodyfat must contain the original complete 71-row dataset,
# with positive DEXfat and the expected numeric predictors.
application <- cis_bodyfat(
  data = original_bodyfat,
  output_dir = file.path(getwd(), "bodyfat_CustomInput")
)

Estimators and calibration

Code label

Manuscript label

Description

MLE

MLE

Unpenalized Gamma regression with log link.

GammaRidge

GR

Mean-weighted quadratic shrinkage of the full-data MLE with HKB tuning.

RefitRidgeHKB

RHKB

HKB tuning applied to the screened reference fit and screened working matrix.

GammaCIS

CIS

Central-scale calibration applied to the same screened reference fit and working matrix.

Refitting is an internal step, not a separate comparison estimator. RHKB is a controlled reference procedure used to isolate the tuning contribution; it is not presented here as an independently published robust estimator. The scripts do not implement the published Mallows-based robust Gamma ridge estimator.

Screening combines squared robust predictor distances with unscaled Gamma deviance residuals from the initial MLE. Predictor distances use a coordinatewise-median location and MASS::cov.rob(method = "mve", nsamp = 1000) covariance, with a chi-square cutoff at probability 0.975 and p degrees of freedom. The absolute unscaled deviance-residual cutoff is 2.5. The two binary retention indicators are multiplied.

If no more than q + 5 observations survive, the code falls back to the full-data MLE and resets the working weights to one. Otherwise it refits the Gamma model on retained rows. Rank and convergence checks can still fail.

For retained rows, define:

$$
D_C=X^\top\operatorname{diag}(w_i\widehat\mu_{i,\mathrm{ref}})X,
\qquad C=\operatorname{median}{\operatorname{diag}(D_C)},
$$

$$
k_C=\frac{C\widehat\phi_{\mathrm{ref}}}
{|\widehat\beta_{\mathrm{ref}}|2^2+10^{-8}},
\qquad
\widehat\beta_C=(D_C+k_C I_q)^{-1}D_C\widehat\beta{\mathrm{ref}}.
$$

Here phi is estimated by Pearson residual sum of squares divided by residual degrees of freedom. RHKB uses q * phi_ref / sum(beta_ref^2) with the same screened inputs. GR uses the analogous full-data quantities and D_mu = X' diag(mu_MLE) X.

The intercept participates in the central diagonal scale, coefficient norm, and ridge penalty. These mean-weighted matrices are working matrices, not the Fisher information for Gamma regression with a log link. The transformation is a quadratic shrinkage step around a fitted reference, not a direct fit of the exact penalized Gamma likelihood.

Numerical implementation and limitations

Shrinkage is computed using scaled augmented QR. Difficult systems trigger MPFR reconstruction from unchanged inputs, with precision escalation and successive-precision agreement and residual checks. The Gamma MLE first uses IRLS; on failure, a QR-coordinate optimization and Newton refinement attempt to solve the same unpenalized likelihood equations. Unresolved failures are retained and logged.

Exact outputs can vary with R/package versions, numerical tolerances, and recovery policies. Reproduction of every value in an earlier archived workbook is not guaranteed by these scripts. Screening can remove uncontaminated observations, including useful leverage points; the initial MLE can be affected by contamination. Neither screening nor shrinkage guarantees lower risk in every design. Common matrix-scale calibration is not general invariance to predictor transformations, and the automatic penalty is not guaranteed to vanish relative to sample size.

Using the core functions with other data

source("cis_core.R")

# predictors: numeric predictor matrix/data frame; response: positive vector.
design <- cis_make_design(predictors)
fit <- cis_fit(design$X, response, screen_seed = 123L)

fit$beta$GammaCIS
fit$errors
fit$meta$k_CIS

# Conversion for this fixed standardization:
if (!is.null(fit$beta$GammaCIS)) {
  beta_original_units <- cis_original_beta(
    fit$beta$GammaCIS, design$center, design$spread
  )
}

cis_fit() requires a finite numeric design with a unit intercept as its first column and a strictly positive finite response. It does not standardize predictors automatically. Individual failed estimates are returned as NULL, with diagnostic messages in errors; inspect these before using coefficients.

Data reference and citation

The application script identifies the original data source as:

Garcia et al. (2005). Obesity Research, 13, 626–634.
https://doi.org/10.1038/oby.2005.67

The dataset is distributed through the R package TH.data. See the package documentation for full attribution and applicable data-use terms:

help("bodyfat", package = "TH.data")
citation("TH.data")

When using this repository, cite the associated CIS manuscript and the data source. No publication DOI or software release DOI is asserted here; add the final manuscript citation when available.

License

No software license is specified by the three supplied scripts. The repository owner should add an appropriate LICENSE file before granting reuse permissions. Third-party data and dependencies retain their own terms.
