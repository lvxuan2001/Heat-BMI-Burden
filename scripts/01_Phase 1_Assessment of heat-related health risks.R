################################################################################
# Phase 1: Assessment of heat-related health risks
#
# Description: We estimate city-specific heat–health associations using a distributed lag non-linear model (DLNM),
#              then pool them via multivariate meta-analysis, and compute:
#              1) Relative risks (RR) at a high-heat quantile centered at the minimum heat exposure
#              2) Heat-attributable numbers via `attrdl`

# Date: 06-02-2026
# R version: 4.5.2
################################################################################

# ==============================================================================
# 0. ENVIRONMENT SETUP
# ==============================================================================

rm(list = ls())
set.seed(12345)

# Load packages
required_packages <- c("dlnm", "mvmeta", "splines", "tsModel", 
                       "dplyr", "haven", "lubridate", "openxlsx")
for (pkg in required_packages) {
  library(pkg, character.only = TRUE)
}
# Version check
if (packageVersion("dlnm") < "2.2.0") {
  stop("dlnm package version >= 2.2.0 is required")
}

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

CONFIG <- list(
  # DLNM parameters
  varfun = "ns",
  varper = c(10, 75, 90),
  min_lag = 0,
  max_lag = 21,
  lag_knots = logknots(21, 3),
  dfseas = 8,
  
  
  # Disease categories
  diseases = c("totalcase", "E", "I", "J", "N", 
               "female", "male", "lower18", "y18to65", "y65plus"),
  
  # Monte Carlo simulation
  
  nsim = 1000,
  
  # Quantile for heat threshold
  heat_threshold_quantile = 0.95,
  
  # File paths
  input_file = "analysis.dta",
  output_rdata = "phase1_results.RData"
)

# ==============================================================================
# 2. DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load dataset
analysis <- read_dta(CONFIG$input_file)
analysis$date <- as.Date(analysis$date)

# Create city list
city_ids <- as.character(unique(analysis$city))
dlist <- lapply(city_ids, function(x) analysis[analysis$city == x, ])
names(dlist) <- city_ids

# City metadata
cities <- analysis %>%
  select(city, cityname, region, region_code, level) %>%
  distinct() %>%
  arrange(city) %>%
  as.data.frame()

cities$city <- as.character(cities$city)
dlist <- dlist[cities$city]

# ==============================================================================
# 3. STAGE 1: CITY-SPECIFIC DLNM MODELS
# ==============================================================================
# Initialize storage
coef_list <- list()
vcov_list <- list()

for (disease in CONFIG$diseases) {
  k <- tolower(disease)
  coef_list[[paste0("coef_", k)]] <- vector("list", length(dlist))
  vcov_list[[paste0("vcov_", k)]] <- vector("list", length(dlist))
}

# Main loop
time_start <- proc.time()[3]

for (i in seq_along(dlist)) {
  if (i %% 20 == 0) cat(sprintf("Processing city %d/%d...\n", i, length(dlist)))
  
  data <- dlist[[i]]
  
  # Define cross-basis
  argvar <- list(
    fun = CONFIG$varfun,
    knots = quantile(data$mean_EHI, CONFIG$varper / 100, na.rm = TRUE),
    Boundary.knots = range(data$mean_EHI, na.rm = TRUE)
  )
  
  cb <- crossbasis(
    data$mean_EHI,
    lag = c(CONFIG$min_lag, CONFIG$max_lag),
    argvar = argvar,
    arglag = list(fun = "ns", knots = CONFIG$lag_knots)
  )
  
  # Centering value (city-specific median)
  cen <- median(data$mean_EHI, na.rm = TRUE)
  
  # Fit models for each disease
  for (disease in CONFIG$diseases) {
    k <- tolower(disease)
    
    # Model formula
    n_years <- length(unique(year(data$date)))
    df_seasonal <- CONFIG$dfseas * n_years
    
    formula_str <- sprintf(
      "%s ~ cb + dow + holiday + PM + O3 + ns(as.numeric(date), df = %d)",
      disease, df_seasonal
    )
    
    model <- tryCatch({
      glm(as.formula(formula_str), data = data, 
          family = quasipoisson, na.action = "na.exclude")
    }, error = function(e) NULL)
    
    if (is.null(model)) {
      warning(sprintf("City %s, disease %s: model fitting failed.", 
                      cities$city[i], disease))
      next
    }
    
    # Reduce to overall cumulative effect
    red <- crossreduce(cb, model, model.link = "log", cen = cen)
    
    coef_list[[paste0("coef_", k)]][[i]] <- coef(red)
    vcov_list[[paste0("vcov_", k)]][[i]] <- vcov(red)
  }
}

time_elapsed <- proc.time()[3] - time_start

# Convert coefficient lists to matrices
for (disease in CONFIG$diseases) {
  k <- tolower(disease)
  coef_list[[paste0("coef_", k)]] <- do.call(rbind, coef_list[[paste0("coef_", k)]])
  rownames(coef_list[[paste0("coef_", k)]]) <- cities$city
}

# ==============================================================================
# 4. STAGE 2: MULTIVARIATE META-ANALYSIS
# ==============================================================================

# Meta-predictors
cities$avgtmean <- sapply(dlist, function(x) mean(x$mean_EHI, na.rm = TRUE))
cities$rangetmean <- sapply(dlist, function(x) diff(range(x$mean_EHI, na.rm = TRUE)))

# Run meta-analysis for each disease
mv_list <- setNames(lapply(CONFIG$diseases, function(disease) {
  k <- tolower(disease)
  
  coef_mat <- coef_list[[paste0("coef_", k)]]
  vcov_lst <- vcov_list[[paste0("vcov_", k)]]
  
  # Remove cities with missing coefficients
  valid_idx <- !sapply(1:nrow(coef_mat), function(i) any(is.na(coef_mat[i, ])))
  mvmeta(
    coef_mat[valid_idx, ] ~ avgtmean + rangetmean,
    S = vcov_lst[valid_idx],
    data = cities[valid_idx, ],
    control = list(showiter = FALSE),
    method = "reml"
  )
}), paste0("mv_", tolower(CONFIG$diseases)))

# Wald test for effect modifiers
fwald <- function(model, var) {
  if (is.null(model)) return(NA)
  ind <- grep(var, names(coef(model)))
  if (length(ind) == 0) return(NA)
  coef_temp <- coef(model)[ind]
  vcov_temp <- vcov(model)[ind, ind]
  waldstat <- coef_temp %*% solve(vcov_temp) %*% coef_temp
  df <- length(coef_temp)
  return(1 - pchisq(waldstat, df))
}

wald_results <- lapply(CONFIG$diseases, function(disease) {
  k <- tolower(disease)
  model <- mv_list[[paste0("mv_", k)]]
  list(
    avgtmean = fwald(model, "avgtmean"),
    rangetmean = fwald(model, "rangetmean")
  )
})
names(wald_results) <- paste0("wald_", tolower(CONFIG$diseases))


# ==============================================================================
# 5. BEST LINEAR UNBIASED PREDICTIONS (BLUPs)
# ==============================================================================

blup_list <- setNames(lapply(CONFIG$diseases, function(disease) {
  k <- tolower(disease)
  model <- mv_list[[paste0("mv_", k)]]
  if (is.null(model)) return(NULL)
  predict(model, vcov = TRUE)
}), paste0("blup_", tolower(CONFIG$diseases)))


# ==============================================================================
# 6. RELATIVE RISK CALCULATION
# ==============================================================================

#' Calculate risk metrics with MHE-centered approach
#' @param coef Coefficient vector from BLUP
#' @param vcov Variance-covariance matrix
#' @param exposure_data Exposure vector for basis construction
#' @param label_info List with Level, Location, Region, Disease
#' @param p_target Target quantile for RR calculation (default: 0.95)
#' @return Data frame with RR estimates and CIs
get_risk_metrics <- function(coef, vcov, exposure_data, label_info, p_target = 0.95) {
  
  # Define basis (consistent with Stage 1)
  argvar <- list(
    fun = CONFIG$varfun,
    knots = quantile(exposure_data, CONFIG$varper / 100, na.rm = TRUE),
    Boundary.knots = range(exposure_data, na.rm = TRUE)
  )
  
  bvar <- onebasis(exposure_data, fun = CONFIG$varfun, 
                   knots = argvar$knots, 
                   Boundary.knots = argvar$Boundary.knots)
  
  # Find minimum heat exposure (MHE) within 1%-99% range
  pr_range <- quantile(exposure_data, c(0.01, 0.99), na.rm = TRUE)
  pred_grid <- seq(pr_range[1], pr_range[2], length.out = 5000)
  
  pred_temp <- crosspred(bvar, coef = coef, vcov = vcov, model.link = "log",
                         at = pred_grid, cen = median(exposure_data, na.rm = TRUE))
  
  mhe_val <- pred_temp$predvar[which.min(pred_temp$allRRfit)]
  
  # Calculate RR at target quantile, centered at MHE
  
  target_val <- quantile(exposure_data, p_target, na.rm = TRUE)
  
  pred_final <- crosspred(bvar, coef = coef, vcov = vcov, model.link = "log",
                          at = target_val, cen = mhe_val)
  
  data.frame(
    Level = label_info$Level,
    Location = label_info$Location,
    Region = label_info$Region,
    Disease = label_info$Disease,
    MHE = mhe_val,
    Target_EHI = as.numeric(target_val),
    RR = as.numeric(pred_final$allRRfit),
    RR_Low = as.numeric(pred_final$allRRlow),
    RR_High = as.numeric(pred_final$allRRhigh),
    stringsAsFactors = FALSE
  )
}

# Calculate RRs
final_results_rr <- data.frame()

# City-level RRs
cat("  Processing city-level RRs...\n")
for (i in seq_along(dlist)) {
  if (i %% 50 == 0) cat(sprintf("    City %d/%d\n", i, length(dlist)))
  
  city_data <- dlist[[i]]$mean_EHI
  
  for (disease in CONFIG$diseases) {
    k <- tolower(disease)
    
    blup_fit <- blup_list[[paste0("blup_", k)]]
    if (is.null(blup_fit)) next
    
    blup_coef <- blup_fit[[i]]$fit
    blup_vcov <- blup_fit[[i]]$vcov
    
    if (is.null(blup_coef)) next
    
    row <- get_risk_metrics(
      coef = blup_coef,
      vcov = blup_vcov,
      exposure_data = city_data,
      label_info = list(
        Level = "City",
        Location = cities$cityname[i],
        Region = cities$region[i],
        Disease = disease
      )
    )
    final_results_rr <- rbind(final_results_rr, row)
  }
}

# National-level RRs
national_EHI <- unlist(lapply(dlist, function(x) x$mean_EHI))

for (disease in CONFIG$diseases) {
  k <- tolower(disease)
  
  model <- mv_list[[paste0("mv_", k)]]
  if (is.null(model)) next
  
  newdata <- data.frame(
    avgtmean = mean(cities$avgtmean),
    rangetmean = mean(cities$rangetmean)
  )
  
  pooled <- predict(model, newdata = newdata, vcov = TRUE)
  
  p_fit <- if (is.list(pooled) && !is.null(pooled$fit)) pooled$fit else pooled[[1]]$fit
  p_vcov <- if (is.list(pooled) && !is.null(pooled$vcov)) pooled$vcov else pooled[[1]]$vcov
  
  row <- get_risk_metrics(
    coef = p_fit,
    vcov = p_vcov,
    exposure_data = national_EHI,
    label_info = list(
      Level = "National",
      Location = "China",
      Region = "National",
      Disease = disease
    )
  )
  final_results_rr <- rbind(final_results_rr, row)
}

# ==============================================================================
# 7. ATTRIBUTABLE FRACTION CALCULATION
# ==============================================================================

source("attrdl.R")

# Initialize storage
tothosp_list <- setNames(lapply(CONFIG$diseases, function(disease) {
  setNames(rep(NA_real_, nrow(cities)), cities$city)
}), paste0("tothosp_", tolower(CONFIG$diseases)))

matsim_list <- setNames(lapply(CONFIG$diseases, function(disease) {
  matrix(NA_real_, nrow(cities), 1, dimnames = list(cities$city, "heat"))
}), paste0("matsim_", tolower(CONFIG$diseases)))

arraysim_list <- setNames(lapply(CONFIG$diseases, function(disease) {
  array(NA_real_, dim = c(nrow(cities), 1, CONFIG$nsim),
        dimnames = list(cities$city, "heat", NULL))
}), paste0("arraysim_", tolower(CONFIG$diseases)))

# Main AF calculation loop
for (i in seq_along(dlist)) {
  if (i %% 20 == 0) cat(sprintf("    City %d/%d\n", i, length(dlist)))
  
  data <- dlist[[i]]
  
  # Define cross-basis (consistent with Stage 1)
  argvar <- list(
    fun = CONFIG$varfun,
    knots = quantile(data$mean_EHI, CONFIG$varper / 100, na.rm = TRUE),
    Boundary.knots = range(data$mean_EHI, na.rm = TRUE)
  )
  
  cb <- crossbasis(
    data$mean_EHI,
    lag = c(CONFIG$min_lag, CONFIG$max_lag),
    argvar = argvar,
    arglag = list(fun = "ns", knots = CONFIG$lag_knots)
  )
  
  # One-basis for MHE calculation
  bvar_mhe <- onebasis(data$mean_EHI, fun = CONFIG$varfun,
                       knots = argvar$knots,
                       Boundary.knots = argvar$Boundary.knots)
  
  for (disease in CONFIG$diseases) {
    k <- tolower(disease)
    
    blup_fit <- blup_list[[paste0("blup_", k)]]
    if (is.null(blup_fit) || is.null(blup_fit[[i]])) next
    
    blup_coef <- blup_fit[[i]]$fit
    blup_vcov <- blup_fit[[i]]$vcov
    
    # Find MHE
    pr <- quantile(data$mean_EHI, c(0.01, 0.99), na.rm = TRUE)
    pred_grid <- seq(pr[1], pr[2], length.out = 5000)
    
    pred_temp <- crosspred(bvar_mhe, coef = blup_coef, vcov = blup_vcov,
                           model.link = "log", at = pred_grid,
                           cen = median(data$mean_EHI, na.rm = TRUE))
    
    mhe_val <- pred_temp$predvar[which.min(pred_temp$allRRfit)]
    
    # Heat threshold (P95)
    ehi_95th <- quantile(data$mean_EHI, CONFIG$heat_threshold_quantile, na.rm = TRUE)
    
    # Point estimate
    matsim_list[[paste0("matsim_", k)]][i, "heat"] <- attrdl(
      data$mean_EHI, cb, data[[disease]],
      coef = blup_coef, vcov = blup_vcov,
      type = "an", dir = "forw",
      cen = mhe_val, range = c(ehi_95th, 100)
    )
    
    # Monte Carlo simulation for CI
    arraysim_list[[paste0("arraysim_", k)]][i, "heat", ] <- attrdl(
      data$mean_EHI, cb, data[[disease]],
      coef = blup_coef, vcov = blup_vcov,
      type = "an", dir = "forw",
      cen = mhe_val, range = c(ehi_95th, 100),
      sim = TRUE, nsim = CONFIG$nsim
    )
    
    # Total hospitalizations
    tothosp_list[[paste0("tothosp_", k)]][i] <- sum(data[[disease]], na.rm = TRUE)
  }
}