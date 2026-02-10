################################################################################
# Phase 2: Health insurance burden assessment
#
# Description: Estimation of heat-attributable inpatient reimbursements (absolute burden) 
#              and their share of total fund expenditure (relative burden) with 
#              Monte Carlo uncertainty propagation, LMDI attribution,variance decomposition,
#              and Compound annual growth rate of GDP growth and absolute burden 
#
# Date: 06-02-2026
# R version: 4.5.2
################################################################################

# ==============================================================================
# 0. ENVIRONMENT SETUP
# ==============================================================================

rm(list = ls())
set.seed(12345)

# Load packages
required_packages <- c("dplyr", "tidyr", "openxlsx", "boot", "haven")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

P2_CONFIG <- list(
  # Monte Carlo parameters
  nsim = 1000,
  conf_level = 0.95,
  zero_floor = 1e-9,
  
  # Analysis periods
  years = c(2020, 2021, 2022, 2023),
  periods = list(c(2020, 2023)),
  
  # Disease categories for burden analysis
  disease_codes = c("I", "E", "J", "N"),
  
  # Bootstrap parameters
  n_boot = 1000,
  
  # File paths
  phase1_input = "phase1_results.RData",
  burden_data_file = "burden_data.dta",
  gdp_pop_file = "GDP_POP.dta"
)

# ==============================================================================
# 2. DATA LOADING
# ==============================================================================
# Load Phase 1 results
if (!file.exists(P2_CONFIG$phase1_input)) {
  stop(sprintf("Phase 1 results file '%s' not found. Run Phase 1 first.", 
               P2_CONFIG$phase1_input))
}
load(P2_CONFIG$phase1_input)


# Load burden data
if (!file.exists(P2_CONFIG$burden_data_file)) {
  stop(sprintf("Burden data file '%s' not found.", P2_CONFIG$burden_data_file))
}
burden_data <- read_dta(P2_CONFIG$burden_data_file)

# Load GDP data (for GDP growth rate analysis)
if (!file.exists(P2_CONFIG$gdp_pop_file)) {
  warning(sprintf("GDP data file '%s' not found. GDP growth analysis will be skipped.", 
                  P2_CONFIG$gdp_pop_file))
  gdp_data <- NULL
} else {
  gdp_data <- read_dta(P2_CONFIG$gdp_pop_file)
  cat(sprintf("  GDP data loaded: %d observations.\n", nrow(gdp_data)))
}

# Validate Phase 1 objects
required_objects <- c("cities", "matsim_list", "arraysim_list", "tothosp_list")
missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects) > 0) {
  stop(sprintf("Missing Phase 1 objects: %s", paste(missing_objects, collapse = ", ")))
}

# Ensure city ID compatibility
cities$city <- as.double(cities$city)

# ==============================================================================
# 3. DATA PREPROCESSING
# ==============================================================================

# Filter cities present in burden data
cities_filtered <- cities %>%
  filter(city %in% burden_data$city) %>%
  arrange(city)

city_indices <- match(cities_filtered$city, cities$city)

# Calculate attributable fractions and standard errors for each disease
calc_af_se <- function(matsim, arraysim, tothosp, idx) {
  af_mean <- (matsim[idx, "heat"] / tothosp[idx]) * 100
  af_upper <- apply(arraysim[idx, , , drop = FALSE], c(1, 2), quantile, 0.975)[, "heat"]
  af_lower <- apply(arraysim[idx, , , drop = FALSE], c(1, 2), quantile, 0.025)[, "heat"]
  af_se <- ((af_upper - af_lower) / tothosp[idx] * 100) / (2 * 1.96)
  
  list(mean = af_mean, se = af_se)
}

# Join AF estimates to burden data
af_data <- data.frame(city = cities_filtered$city)

for (code in P2_CONFIG$disease_codes) {
  k <- tolower(code)
  
  if (!paste0("matsim_", k) %in% names(matsim_list)) {
    warning(sprintf("Disease %s not found in Phase 1 results, skipping.", code))
    next
  }
  
  af_result <- calc_af_se(
    matsim_list[[paste0("matsim_", k)]],
    arraysim_list[[paste0("arraysim_", k)]],
    tothosp_list[[paste0("tothosp_", k)]],
    city_indices
  )
  
  af_data[[paste0(code, "Heat")]] <- af_result$mean
  af_data[[paste0(code, "Heat_SE")]] <- af_result$se
}

burden_data <- burden_data %>%
  left_join(af_data, by = "city")

# Convert fund units 
burden_data <- burden_data %>%
  mutate(
    Fund_Total_Yuan = Healthcare_Fund * 10000,
    Fund_Emp_Yuan = Fund_Emp * 10000,
    Fund_Res_Yuan = Fund_Res * 10000,
    cityname = as.character(cityname),
    Region = as.character(Region)
  )

# Merge GDP data if available
if (!is.null(gdp_data)) {
  gdp_data <- gdp_data %>%
    mutate(city = as.double(city)) %>%
    select(city, year, any_of(c("GDP", "gdp", "GDP_pc", "gdp_pc")))
  
  # Standardize GDP column name
  if ("gdp" %in% names(gdp_data) && !"GDP" %in% names(gdp_data)) {
    gdp_data <- gdp_data %>% rename(GDP = gdp)
  }
  if ("gdp_pc" %in% names(gdp_data) && !"GDP_pc" %in% names(gdp_data)) {
    gdp_data <- gdp_data %>% rename(GDP_pc = gdp_pc)
  }
  
  burden_data <- burden_data %>%
    left_join(gdp_data, by = c("city", "year"))
}

# ==============================================================================
# 4. SIMULATION ENGINE
# ==============================================================================

n_rows <- nrow(burden_data)

#' Generate simulation matrix with truncation at zero
#' @param mean_vec Vector of means
#' @param se_vec Vector of standard errors
#' @param n_sim Number of simulations
#' @return Matrix (n_obs x n_sim) of simulated values
gen_sim_matrix <- function(mean_vec, se_vec, n_sim = P2_CONFIG$nsim) {
  mat <- matrix(rnorm(length(mean_vec) * n_sim, mean = mean_vec, sd = se_vec), 
                nrow = length(mean_vec))
  mat[mat < 0] <- 0
  mat
}

# Generate AF simulation matrices
mat_AF <- list()
for (code in P2_CONFIG$disease_codes) {
  col_mean <- paste0(code, "Heat")
  col_se <- paste0(code, "Heat_SE")
  
  if (!col_mean %in% names(burden_data)) next
  
  mat_AF[[code]] <- gen_sim_matrix(
    burden_data[[col_mean]], 
    burden_data[[col_se]]
  )
}

#' Calculate expenditure matrix
#' @param AF_mat Attributable fraction matrix
#' @param R_vec Hospitalization rate
#' @param Exp_vec Per-capita expenditure
#' @param Ins_vec Insurance coverage rate
#' @param Rate_vec Reimbursement rate
#' @param Pop_vec Population
#' @return Expenditure matrix
calc_exp_matrix <- function(AF_mat, R_vec, Exp_vec, Ins_vec, Rate_vec, Pop_vec) {
  coef <- (R_vec * Pop_vec * Ins_vec * Exp_vec * Rate_vec) / 100
  AF_mat * coef
}

# Calculate expenditure matrices by insurance type and disease
Exp_matrices <- list()

for (ins_type in c("Emp", "Res")) {
  ins_rate_col <- paste0("InsuranceRate_", tolower(ins_type))
  rate_col <- paste0(tolower(ins_type), "_rate")
  
  for (code in P2_CONFIG$disease_codes) {
    if (!code %in% names(mat_AF)) next
    
    key <- paste0(code, "_", ins_type)
    Exp_matrices[[key]] <- calc_exp_matrix(
      mat_AF[[code]],
      burden_data[[paste0("R_", code)]],
      burden_data[[paste0("Expenditure_", code)]],
      burden_data[[ins_rate_col]],
      burden_data[[rate_col]],
      burden_data$pop
    )
  }
}

# Aggregate by insurance type
Mat_Total_Emp <- Reduce(`+`, Exp_matrices[grep("_Emp$", names(Exp_matrices))])
Mat_Total_Res <- Reduce(`+`, Exp_matrices[grep("_Res$", names(Exp_matrices))])
Mat_Total_All <- Mat_Total_Emp + Mat_Total_Res

# Create disease-specific totals
Exp_All <- list()
for (code in P2_CONFIG$disease_codes) {
  emp_key <- paste0(code, "_Emp")
  res_key <- paste0(code, "_Res")
  if (emp_key %in% names(Exp_matrices) && res_key %in% names(Exp_matrices)) {
    Exp_All[[code]] <- Exp_matrices[[emp_key]] + Exp_matrices[[res_key]]
  }
}

# Clean up AF matrices
rm(mat_AF)
gc()

# ==============================================================================
# 5. BURDEN ANALYSIS
# ==============================================================================

cat("\n========== Task 1: Burden Analysis ==========\n")

#' Aggregate expenditures and calculate statistics
#' @param input_mat Expenditure matrix (n_obs x n_sim)
#' @param group_vec Grouping vector
#' @param fund_vec Fund values for relative burden
#' @return Data frame with absolute and relative burden statistics
aggregate_and_stats <- function(input_mat, group_vec, fund_vec) {
  agg_exp <- rowsum(input_mat, group = group_vec, reorder = TRUE)
  agg_fund <- tapply(fund_vec, group_vec, sum, na.rm = TRUE)
  agg_fund <- agg_fund[rownames(agg_exp)]
  
  # Absolute burden (million yuan)
  exp_mean <- rowMeans(agg_exp) / 1e6
  exp_low <- apply(agg_exp, 1, quantile, 0.025) / 1e6
  exp_high <- apply(agg_exp, 1, quantile, 0.975) / 1e6
  
  # Relative burden (per mille of fund)
  rel_mat <- agg_exp / as.vector(agg_fund) * 1000
  rel_mean <- rowMeans(rel_mat)
  rel_low <- apply(rel_mat, 1, quantile, 0.025)
  rel_high <- apply(rel_mat, 1, quantile, 0.975)
  
  data.frame(
    Group = rownames(agg_exp),
    Abs_Mean = exp_mean,
    Abs_Low = exp_low,
    Abs_High = exp_high,
    Rel_Mean = rel_mean,
    Rel_Low = rel_low,
    Rel_High = rel_high,
    stringsAsFactors = FALSE
  )
}

# Define scenarios and dimensions
scenarios <- list(
  Total = list(
    fund = burden_data$Fund_Total_Yuan,
    mats = c(list(Total = Mat_Total_All), Exp_All)
  ),
  Employee = list(
    fund = burden_data$Fund_Emp_Yuan,
    mats = list(Total = Mat_Total_Emp)
  ),
  Resident = list(
    fund = burden_data$Fund_Res_Yuan,
    mats = list(Total = Mat_Total_Res)
  )
)

dims <- list(
  City = burden_data$cityname,
  Level = burden_data$level,
  National = rep("China", n_rows)
)

# Run burden analysis
results_burden <- data.frame()

for (ins_type in names(scenarios)) {
  cat(sprintf("  Processing: %s...\n", ins_type))
  
  curr_fund <- scenarios[[ins_type]]$fund
  disease_mats <- scenarios[[ins_type]]$mats
  
  for (disease_name in names(disease_mats)) {
    curr_mat <- disease_mats[[disease_name]]
    if (is.null(curr_mat)) next
    
    for (dim_name in names(dims)) {
      curr_group <- dims[[dim_name]]
      
      # By year
      for (y in P2_CONFIG$years) {
        idx <- burden_data$year == y
        if (sum(idx) == 0) next
        
        res <- aggregate_and_stats(curr_mat[idx, , drop = FALSE], 
                                   curr_group[idx], curr_fund[idx])
        res$Dimension <- dim_name
        res$Year <- as.character(y)
        res$Insurance <- ins_type
        res$Disease <- disease_name
        results_burden <- rbind(results_burden, res)
      }
      
      # Total period
      res_total <- aggregate_and_stats(curr_mat, curr_group, curr_fund)
      res_total$Dimension <- dim_name
      res_total$Year <- "2020-2023"
      res_total$Insurance <- ins_type
      res_total$Disease <- disease_name
      results_burden <- rbind(results_burden, res_total)
    }
  }
}

# ==============================================================================
# 6. LMDI ATTRIBUTION ANALYSIS
# ==============================================================================

#' Logarithmic mean
L_mean <- function(x, y) {
  res <- (x - y) / (log(x) - log(y))
  res[x == y] <- x[x == y]
  res[!is.finite(res)] <- 0
  res
}

#' LMDI attribution decomposition
#' @param mat_emp_s Employee expenditure matrix (start period)
#' @param mat_emp_e Employee expenditure matrix (end period)
#' @param fund_emp_s Employee fund (start)
#' @param fund_emp_e Employee fund (end)
#' @param mat_res_s Resident expenditure matrix (start)
#' @param mat_res_e Resident expenditure matrix (end)
#' @param fund_res_s Resident fund (start)
#' @param fund_res_e Resident fund (end)
#' @return Data frame with attribution results
calc_lmdi_attribution <- function(mat_emp_s, mat_emp_e, fund_emp_s, fund_emp_e,
                                  mat_res_s, mat_res_e, fund_res_s, fund_res_e) {
  
  E_emp_0 <- colSums(mat_emp_s, na.rm = TRUE)
  E_emp_t <- colSums(mat_emp_e, na.rm = TRUE)
  F_emp_0 <- sum(fund_emp_s, na.rm = TRUE)
  F_emp_t <- sum(fund_emp_e, na.rm = TRUE)
  
  E_res_0 <- colSums(mat_res_s, na.rm = TRUE)
  E_res_t <- colSums(mat_res_e, na.rm = TRUE)
  F_res_0 <- sum(fund_res_s, na.rm = TRUE)
  F_res_t <- sum(fund_res_e, na.rm = TRUE)
  
  E_tot_0 <- E_emp_0 + E_res_0
  E_tot_t <- E_emp_t + E_res_t
  F_tot_0 <- max(F_emp_0 + F_res_0, P2_CONFIG$zero_floor)
  F_tot_t <- max(F_emp_t + F_res_t, P2_CONFIG$zero_floor)
  
  R_tot_0 <- E_tot_0 / F_tot_0
  R_tot_t <- E_tot_t / F_tot_t
  W_tot <- L_mean(R_tot_t, R_tot_0)
  
  L_E_tot <- L_mean(E_tot_t, E_tot_0)
  L_E_tot[L_E_tot == 0] <- P2_CONFIG$zero_floor
  
  # Cost effects
  Eff_Cost_Emp <- W_tot * ((E_emp_t - E_emp_0) / L_E_tot)
  Eff_Cost_Res <- W_tot * ((E_res_t - E_res_0) / L_E_tot)
  
  # Fund effects
  L_F_tot <- max(L_mean(F_tot_t, F_tot_0), P2_CONFIG$zero_floor)
  Eff_Fund_Emp <- W_tot * ((F_emp_0 - F_emp_t) / L_F_tot)
  Eff_Fund_Res <- W_tot * ((F_res_0 - F_res_t) / L_F_tot)
  
  Total_Change <- R_tot_t - R_tot_0
  
  # Summary statistics (convert to per mille)
  get_stats <- function(vec, prefix) {
    vec <- vec * 1000
    data.frame(
      Mean = mean(vec, na.rm = TRUE),
      Low = quantile(vec, 0.025, na.rm = TRUE),
      High = quantile(vec, 0.975, na.rm = TRUE)
    ) %>% setNames(paste0(prefix, c("_Mean", "_Low", "_High")))
  }
  
  cbind(
    get_stats(Eff_Cost_Emp, "Emp_Cost"),
    get_stats(Eff_Fund_Emp, "Emp_Fund"),
    get_stats(Eff_Cost_Res, "Res_Cost"),
    get_stats(Eff_Fund_Res, "Res_Fund"),
    get_stats(Total_Change, "Total")
  )
}

# Run LMDI analysis
results_attribution <- data.frame()

for (dim_name in names(dims)) {
  unique_groups <- unique(dims[[dim_name]])
  
  for (grp in unique_groups) {
    for (p in P2_CONFIG$periods) {
      y_start <- p[1]
      y_end <- p[2]
      
      idx_s <- (dims[[dim_name]] == grp) & (burden_data$year == y_start)
      idx_e <- (dims[[dim_name]] == grp) & (burden_data$year == y_end)
      
      if (sum(idx_s) == 0 || sum(idx_e) == 0) next
      
      res <- calc_lmdi_attribution(
        Mat_Total_Emp[idx_s, , drop = FALSE], Mat_Total_Emp[idx_e, , drop = FALSE],
        burden_data$Fund_Emp_Yuan[idx_s], burden_data$Fund_Emp_Yuan[idx_e],
        Mat_Total_Res[idx_s, , drop = FALSE], Mat_Total_Res[idx_e, , drop = FALSE],
        burden_data$Fund_Res_Yuan[idx_s], burden_data$Fund_Res_Yuan[idx_e]
      )
      
      res$Dimension <- dim_name
      res$Group <- grp
      res$Period <- paste0(y_start, "-", y_end)
      
      results_attribution <- rbind(results_attribution, res)
    }
  }
}

# ==============================================================================
# 7. BURDEN CAGR ANALYSIS 
# ==============================================================================

#' Calculate CAGR statistics (Monte Carlo)
#' Formula: (End/Start)^(1/n) - 1
calc_cagr_stats <- function(mat_end, mat_start, n_years) {
  # Sum columns to get total cost per simulation
  total_end <- colSums(mat_end, na.rm = TRUE)
  total_start <- colSums(mat_start, na.rm = TRUE)
  
  # Avoid division by zero or negative base
  total_start[total_start <= 1e-6] <- NA
  
  # Calculate CAGR vector
  cagr_vec <- (total_end / total_start)^(1 / n_years) - 1
  cagr_vec <- cagr_vec * 100 # Convert to percentage
  
  data.frame(
    CAGR_Mean = mean(cagr_vec, na.rm = TRUE),
    CAGR_Low = quantile(cagr_vec, 0.025, na.rm = TRUE),
    CAGR_High = quantile(cagr_vec, 0.975, na.rm = TRUE)
  )
}

results_growth <- data.frame()
years_sorted <- sort(P2_CONFIG$years)
start_year <- min(years_sorted) # 2020
end_year <- max(years_sorted)   # 2023
n_years <- end_year - start_year # 3

for (ins_type in names(scenarios)) {
  disease_mats <- scenarios[[ins_type]]$mats
  
  for (disease_name in names(disease_mats)) {
    curr_mat <- disease_mats[[disease_name]]
    if (is.null(curr_mat)) next
    
    for (dim_name in names(dims)) {
      curr_group <- dims[[dim_name]]
      unique_groups <- unique(curr_group)
      
      for (grp in unique_groups) {
        # Only calculate Overall Period (2023 vs 2020)
        idx_end <- (curr_group == grp) & (burden_data$year == end_year)
        idx_start <- (curr_group == grp) & (burden_data$year == start_year)
        
        if (sum(idx_end) > 0 && sum(idx_start) > 0) {
          stats <- calc_cagr_stats(curr_mat[idx_end, , drop = FALSE], 
                                   curr_mat[idx_start, , drop = FALSE],
                                   n_years)
          stats$Dimension <- dim_name
          stats$Group <- grp
          stats$Insurance <- ins_type
          stats$Disease <- disease_name
          stats$Period <- paste0(start_year, "-", end_year)
          
          results_growth <- rbind(results_growth, stats)
        }
      }
    }
  }
}

# ==============================================================================
# 8. GDP CAGR ANALYSIS
# ==============================================================================

results_gdp_growth <- data.frame()

if (!is.null(gdp_data) && "GDP" %in% names(burden_data)) {
  
  calc_gdp_cagr <- function(gdp_end, gdp_start, n_years) {
    if (gdp_start <= 1e-6) return(NA)
    ((gdp_end / gdp_start)^(1 / n_years) - 1) * 100
  }
  
  for (dim_name in names(dims)) {
    curr_group <- dims[[dim_name]]
    unique_groups <- unique(curr_group)
    
    for (grp in unique_groups) {
      # Only calculate Overall Period (2023 vs 2020)
      idx_end <- (curr_group == grp) & (burden_data$year == end_year)
      idx_start <- (curr_group == grp) & (burden_data$year == start_year)
      
      if (sum(idx_end) > 0 && sum(idx_start) > 0) {
        gdp_end <- sum(burden_data$GDP[idx_end], na.rm = TRUE)
        gdp_start <- sum(burden_data$GDP[idx_start], na.rm = TRUE)
        
        if (gdp_end > 0 && gdp_start > 0) {
          cagr <- calc_gdp_cagr(gdp_end, gdp_start, n_years)
          
          results_gdp_growth <- rbind(results_gdp_growth, data.frame(
            Dimension = dim_name,
            Group = grp,
            Period = paste0(start_year, "-", end_year),
            GDP_CAGR = cagr,
            GDP_Curr = gdp_end,
            GDP_Prev = gdp_start,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  cat(sprintf("  GDP CAGR analysis complete: %d records.\n", nrow(results_gdp_growth)))
  
} else {
  cat("  GDP data not available, skipping GDP analysis.\n")
}
# ==============================================================================
# 9. FIELDS VARIANCE DECOMPOSITION
# ==============================================================================

# Prepare decomposition variables
decomp_data <- burden_data %>%
  mutate(
    AB_Total_Mean = rowMeans(Mat_Total_All),
    RB_Total = (AB_Total_Mean / Fund_Total_Yuan) * 1000,
    
    R_sum = R_I + R_E + R_J + R_N,
    R_sum_safe = pmax(R_sum, P2_CONFIG$zero_floor),
    
    AF_weighted = (IHeat * R_I + EHeat * R_E + JHeat * R_J + NHeat * R_N) / R_sum_safe,
    Exp_weighted = (Expenditure_I * R_I + Expenditure_E * R_E + 
                      Expenditure_J * R_J + Expenditure_N * R_N) / R_sum_safe,
    
    IR_total = InsuranceRate_emp + InsuranceRate_res,
    IR_total_safe = pmax(IR_total, P2_CONFIG$zero_floor),
    theta_weighted = (emp_rate * InsuranceRate_emp + res_rate * InsuranceRate_res) / IR_total_safe,
    
    Fund = Fund_Total_Yuan,
    POP = pop,
    R_total = R_sum,
    
    # Log transformations
    ln_RB = log(pmax(RB_Total, P2_CONFIG$zero_floor)),
    ln_AF = log(pmax(AF_weighted, P2_CONFIG$zero_floor)),
    ln_R = log(pmax(R_total, P2_CONFIG$zero_floor)),
    ln_POP = log(pmax(POP, P2_CONFIG$zero_floor)),
    ln_IR = log(pmax(IR_total, P2_CONFIG$zero_floor)),
    ln_Exp = log(pmax(Exp_weighted, P2_CONFIG$zero_floor)),
    ln_theta = log(pmax(theta_weighted, P2_CONFIG$zero_floor)),
    ln_Fund = log(pmax(Fund, P2_CONFIG$zero_floor))
  ) %>%
  filter(
    is.finite(ln_RB), is.finite(ln_AF), is.finite(ln_Fund),
    RB_Total > 0, AF_weighted > 0
  )

#' Fields variance decomposition
#' @param data Data frame with log-transformed variables
#' @param year_filter Optional year filter
#' @param group_filter Optional list of group filters
#' @param n_boot Number of bootstrap iterations
#' @return Data frame with decomposition results
fields_decomposition <- function(data, year_filter = NULL, group_filter = NULL,
                                 n_boot = P2_CONFIG$n_boot) {
  
  if (!is.null(year_filter)) data <- data %>% filter(year %in% year_filter)
  if (!is.null(group_filter)) {
    for (col in names(group_filter)) {
      if (col %in% names(data)) {
        data <- data %>% filter(.data[[col]] == group_filter[[col]])
      }
    }
  }
  
  if (nrow(data) < 10) return(NULL)
  
  components <- data %>%
    transmute(
      Y = ln_RB,
      X_AF = ln_AF,
      X_R = ln_R,
      X_POP = ln_POP,
      X_IR = ln_IR,
      X_Exp = ln_Exp,
      X_theta = ln_theta,
      X_Fund = -ln_Fund
    ) %>%
    filter(complete.cases(.)) %>%
    as.data.frame()
  
  if (nrow(components) < 10) return(NULL)
  
  # Decomposition function
  calc_decomp <- function(comp_data, indices = NULL) {
    if (!is.null(indices)) comp_data <- comp_data[indices, , drop = FALSE]
    if (nrow(comp_data) < 5) return(rep(NA, 7))
    
    Y <- comp_data$Y
    var_Y <- var(Y, na.rm = TRUE)
    if (is.na(var_Y) || var_Y < 1e-10) return(rep(NA, 7))
    
    factors <- c("X_AF", "X_R", "X_POP", "X_IR", "X_Exp", "X_theta", "X_Fund")
    sapply(factors, function(f) {
      cov_val <- cov(comp_data[[f]], Y, use = "complete.obs")
      if (is.na(cov_val)) NA else cov_val / var_Y
    })
  }
  
  point_est <- calc_decomp(components)
  if (all(is.na(point_est))) return(NULL)
  
  # Bootstrap CI
  boot_results <- tryCatch({
    boot(data = components, statistic = function(d, i) calc_decomp(d, i), R = n_boot)
  }, error = function(e) NULL)
  
  alpha <- 1 - P2_CONFIG$conf_level
  if (!is.null(boot_results)) {
    ci_low <- apply(boot_results$t, 2, quantile, probs = alpha/2, na.rm = TRUE)
    ci_high <- apply(boot_results$t, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  } else {
    ci_low <- ci_high <- rep(NA, 7)
  }
  
  data.frame(
    Factor = c("Attributable Fraction", "Hospitalization Rate", "Population",
               "Insurance Coverage", "Per-capita Cost", "Reimbursement Rate", "Fund Scale"),
    Code = c("AF", "R", "POP", "IR", "Exp", "theta", "Fund"),
    Contribution = point_est * 100,
    CI_Low = ci_low * 100,
    CI_High = ci_high * 100,
    Variance_Y = var(components$Y, na.rm = TRUE),
    N = nrow(components),
    stringsAsFactors = FALSE
  )
}

# Run Fields decomposition
results_fields <- data.frame()

# By year
for (y in P2_CONFIG$years) {
  res <- fields_decomposition(decomp_data, year_filter = y)
  if (!is.null(res)) {
    res$Year <- as.character(y)
    res$Dimension <- "National"
    res$Group <- "China"
    results_fields <- rbind(results_fields, res)
  }
}

# Pooled
res_pooled <- fields_decomposition(decomp_data)
if (!is.null(res_pooled)) {
  res_pooled$Year <- "2020-2023"
  res_pooled$Dimension <- "National"
  res_pooled$Group <- "China"
  results_fields <- rbind(results_fields, res_pooled)
}

# By city level
levels_available <- unique(decomp_data$level)
levels_available <- levels_available[!is.na(levels_available)]

for (lv in levels_available) {
  res <- fields_decomposition(decomp_data, group_filter = list(level = lv))
  if (!is.null(res)) {
    res$Year <- "2020-2023"
    res$Dimension <- "City_Level"
    res$Group <- as.character(lv)
    results_fields <- rbind(results_fields, res)
  }
}