# Code to analyze the smoking data from GSE148375 using debiased lasso with multiple splitting.
# Only considers a single set of target covariates for illustration purposes.
# To analyze all SNPs, re-run after setting target_snps to each individual SNP in turn.
# - Note: this can take a while even with parallel processing


library(tidyverse)
library(glmnet)
library(GEOquery)
library(foreach)

# Function to extract patient data from GSM object
get_patient_data <- function(gsm){
  metadata <- Meta(gsm)
  characteristics <- metadata[["characteristics_ch1"]]
  smoking_status <- sub("smoking_status: ", "", characteristics[startsWith(characteristics, "smoking_status: ")])
  age <- sub("age: ", "", characteristics[startsWith(characteristics, "age: ")])
  gender <- sub("gender: ", "", characteristics[startsWith(characteristics, "gender: ")])
  ethnicity <- sub("ethnicity: ", "", characteristics[startsWith(characteristics, "ethnicity: ")])
  
  patient_data <- c(patient_id = metadata[["description"]],
                             smoking_status = ifelse(length(smoking_status) == 0, NA, smoking_status),
                             age = ifelse(length(age) == 0, NA, as.numeric(age)),
                             gender = ifelse(length(gender) == 0, NA, gender),
                             ethnicity = ifelse(length(ethnicity) == 0, NA, ethnicity))
  return(patient_data)
}

# Function to replace genotypes with minor allele counts
minor_allele_count <- function(snp_row){
  genotypes <- unlist(snp_row)
  success_indx <- genotypes != "NC"
  minor_allele_counts <- vector(length = length(genotypes))
  minor_allele_counts[!success_indx] <- NA
  
  # identify alleles
  alleles <- unique(genotypes[success_indx]) %>% 
    paste(sep = "", collapse = "") %>% 
    strsplit(split = "") %>% 
    unlist() %>% 
    unique()
  
  if (length(alleles) == 0){
    if (sum(success_indx) != 0){
      stop("Alleles not identified")
    }
  } 
  else{
    if (length(alleles) == 1){
      minor_allele_counts[success_indx] <- rep(0, sum(success_indx))
    }
    else{
      if (length(alleles) == 2){
        characters <- paste(genotypes[success_indx], sep = "", collapse = "")
        allele_counts <- c(str_count(characters, alleles[1]),
                           str_count(characters, alleles[2]))
        minor_allele <- alleles[allele_counts == min(allele_counts)]
        
        # count occurrence of each allele
        minor_allele_counts[success_indx] <- str_count(genotypes[success_indx], minor_allele[1])
      }
      else{
        stop("More than two alleles identified")
      }
    }
  }
  
  return(minor_allele_counts)
}

# Function to check for HW equilibrium
check_hw_equilibrium <- function(minor_allele_counts, test = "chisq"){
  minor_allele_counts <- unlist(minor_allele_counts)
  success_indx <- !is.na(minor_allele_counts)
  minor_allele_count_frequencies <- c(sum(minor_allele_counts[success_indx] == 0),
                                      sum(minor_allele_counts[success_indx] == 1),
                                      sum(minor_allele_counts[success_indx] == 2))
  p_minor_allele <- (minor_allele_count_frequencies[2] + 2 * minor_allele_count_frequencies[3]) / (2 * sum(minor_allele_count_frequencies))
  
  null_dist <- c((1 - p_minor_allele)^2, 
                 2 * p_minor_allele * (1 - p_minor_allele), 
                 p_minor_allele^2)

  hw_equilibrium_pval <- tryCatch(if (test == "chisq"){
    chisq.test(minor_allele_count_frequencies, p = null_dist)$p.value
  } else{
    #fisher.test(matrix(minor_allele_count_frequencies[c(1, 2, 2, 3)] * c(1, 0.5, 0.5, 1), nrow = 2), 
    #            simulate.p.value = FALSE)$p.value
  },
  error = function(e) NA,
  warning = function(e) NA)
  return(hw_equilibrium_pval)
}

# Load Data
nicotine_data <- read_delim("GSE148375_processed_data.txt")
other_data <- getGEO(filename = "GSE148375_family (2).soft.gz")
snp_info <- Table(GPLList(other_data)[[1]])


# Data Cleaning -----------------------------------------------------------

# 1. exclude insertions/deletions
snps_without_metadata <- nicotine_data[!(nicotine_data$ID_REF %in% snp_info$ID), ]
snp_info <- snp_info[(snp_info$ID %in% nicotine_data$ID_REF) &
                       !(snp_info$SNP %in% c("[I/D]", "[D/I]")), ]

nicotine_data <- nicotine_data[nicotine_data$ID_REF %in% snp_info$ID, ]

# extract non-SNP data for patients
patient_info <- as.data.frame(t(sapply(GSMList(other_data), FUN = get_patient_data)))
patient_info$age <- as.numeric(patient_info$age)

# replace genotypes with number of copies of minor allele
nicotine_allele_counts <- t(apply(nicotine_data[, -1], MARGIN = 1, FUN = minor_allele_count))

# 2. exclude call rate < 95%
high_call_rate_indx <- rowMeans(is.na(nicotine_allele_counts)) < 0.05

# 3. exclude if not in HW equilibrium (chi^2 test p<10^-6 threshold)
hw_equil_pvals <- apply(nicotine_allele_counts, MARGIN = 1, FUN = check_hw_equilibrium)
hw_equil_indx <- hw_equil_pvals >= 10^-6

# 4. remove if minor allele freq < 0.01
overall_minor_allele_freqs <- rowSums(nicotine_allele_counts, na.rm = TRUE) / 
  (2 * rowSums(!is.na(nicotine_allele_counts)))

snp_missing_vals <- nicotine_allele_counts %>% 
  is.na() %>% 
  rowSums()

snp_subset_indx <- high_call_rate_indx & !is.na(hw_equil_indx) & hw_equil_indx & 
  (overall_minor_allele_freqs >= 0.01)
snp_subset_ids <- nicotine_data$ID_REF[snp_subset_indx]

# 5. exclude patients missing more than 1% of remaining SNPs
patient_missing_vals <- nicotine_allele_counts[snp_subset_indx, ] %>% 
  is.na() %>% 
  colSums()

patient_subset_indx <- patient_missing_vals < 0.01 * sum(snp_subset_indx)
patient_subset_ids <- colnames(nicotine_data[, -1])[patient_subset_indx]


# Impute SNPs -------------------------------------------------------------
data_matrix <- t(nicotine_allele_counts[snp_subset_indx, patient_subset_indx])
colnames(data_matrix) <- snp_subset_ids
patient_subset_data <- patient_info[match(patient_subset_ids, patient_info$patient_id), 
                                    c("patient_id", "smoking_status", "age", "gender")]
valid_smoking_indx <- patient_subset_data$smoking_status %in% c("Non-smoker", "Smoker")
data_matrix <- data_matrix[valid_smoking_indx, ]
patient_subset_data <- patient_subset_data[valid_smoking_indx, ]

# for each SNP, impute missing minor allele frequencies from observed data
for (i in seq_len(ncol(data_matrix))){
  missing_data_indx <- is.na(data_matrix[, i])
  data_matrix[missing_data_indx, i] <- sample(data_matrix[!missing_data_indx, i], 
                                              sum(missing_data_indx), 
                                              replace = TRUE)
}

data_matrix <- cbind(patient_id = as.numeric(patient_subset_data$patient_id),
                     smoking_status = ifelse(patient_subset_data$smoking_status == "Smoker", 1, 0),
                     age = patient_subset_data$age,
                     male = ifelse(patient_subset_data$gender == "Male", 1, 0),
                     data_matrix)


# Data analysis -----------------------------------------------------------

# Function to return estimated coefficients from lasso, lambda chosen by cross validation
multisplit_selection <- function(X, y){
  n <- nrow(X)
  subsample_indx <- rep(FALSE, n)
  subsample_indx[sample(1:n, n / 2, replace = FALSE)] <- TRUE
  selection_lasso_fit <- cv.glmnet(X[subsample_indx, ], y[subsample_indx],
                                   family = "binomial", standardize = FALSE)
  lasso_coef <- coef(selection_lasso_fit, s = "lambda.min")
  coef_names <- rownames(lasso_coef)
  lasso_coef <- as.data.frame(t(as.vector(lasso_coef)))
  colnames(lasso_coef) <- coef_names
  results <- list(subsample_indx = as.data.frame(t(subsample_indx)),
                  beta_lasso = lasso_coef)
  return(results)
}

# Standardize age
data_matrix[, "age"] <- (data_matrix[, "age"] - mean(data_matrix[, "age"])) / sd(data_matrix[, "age"])


# Sample splitting with lasso selection -----------------------------------
# For each split, save the sampling indicators and selected model
n_cores <- 10 # number of cores for parallel processing
clust <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clust)

B <- 1000
selection_results <- foreach(i = seq_len(B), .packages = c("glmnet")) %dopar% {
  multisplit_selection(X = data_matrix[, -c(1:2)], y = data_matrix[, "smoking_status"])
}

# Aggregate list of results
all_selection_results <- selection_results[[1]]
for (j in seq_along(all_selection_results)){
  if ("data.frame" %in% class(all_selection_results[[j]])){
    tmp <- as.data.frame(matrix(nrow = length(selection_results) - 1, 
                                ncol = ncol(all_selection_results[[j]])))
    colnames(tmp) <- colnames(all_selection_results[[j]])
    all_selection_results[[j]] <- rbind(all_selection_results[[j]], tmp)
    
  }
}

for (i in seq_along(selection_results)){
  if (i > 1){
    for (j in seq_along(all_selection_results)){
      all_selection_results[[j]][i, ] <- c(selection_results[[i]][[j]])
    }
  }
}

all_selection_results$selected_covariates_indx <- all_selection_results$beta_lasso != 0


# Debiased lasso estimates using selected models --------------------------
target_snps <- snp_info %>% 
  filter(((Chr == 15) & (MapInfo %in% 78857862:78933552)) | # CHRNA5/A3/B4
         ((Chr == 8) & (MapInfo %in% 42552509:42623929)) | # CHRNB3/A6
         ((Chr == 19) & (MapInfo %in% 41171810:41524303)),
         ID %in% colnames(data_matrix))
target_covariates <- c("age", "male", target_snps$ID)[-17] # exclude exact duplicate
target_covariate_indx <- match(target_covariates, colnames(data_matrix))
multisplit_coefficients <- as.data.frame(matrix(nrow = B, ncol = length(target_covariates) + 1))
for (b in seq_len(B)){
  subsample_model_indx <- which(all_selection_results$selected_covariates_indx[b, -1]) + 2
  subsample_model_indx <- unique(c(target_covariate_indx, subsample_model_indx))
  subsample_data_indx <- all_selection_results$subsample_indx[b, ]
  
  # Fit debiased lasso
  subsample_lasso_fit <- cv.glmnet(data_matrix[!subsample_data_indx, subsample_model_indx], 
                                   data_matrix[!subsample_data_indx, "smoking_status"],
                                   family = "binomial", standardize = FALSE)
  subsample_lasso_coef <- as.vector(coef(subsample_lasso_fit, s = "lambda.min"))
  multisplit_coefficients[b, ] <- coef(glm(data_matrix[!subsample_data_indx, "smoking_status"] ~ 
                                             data_matrix[!subsample_data_indx, subsample_model_indx],
                                           family = "binomial", control = list(maxit = 1), 
                                           start = subsample_lasso_coef))[1:(length(target_covariates) + 1)]
}
colnames(multisplit_coefficients) <- c("(Intercept)", target_covariates)


# Multiple splitting estimation results -----------------------------------
coef_ests <- colMeans(multisplit_coefficients)
coef_stderrs <- vector(length = length(coef_ests))
centered_sample_indicators <- scale(!all_selection_results$subsample_indx, center = TRUE, scale = FALSE)
n <- nrow(data_matrix)
n2 <- ceiling(n/2)
for (j in seq_along(coef_ests)){
  cov_hat_j <- t(multisplit_coefficients[, j] - coef_ests[j]) %*% centered_sample_indicators / B 
  V_hat_j <- (n * (n - 1)) * sum(cov_hat_j^2) / (n - n2)^2
  coef_stderrs[j] <- sqrt(max(V_hat_j - ((n / B^2) * (n2 / (n - n2)) * 
                                           sum((multisplit_coefficients[, j] - coef_ests[j])^2)),
                              0))
}

results_table <- data.frame(id = names(coef_ests),
                            coef_est = coef_ests,
                            stderr = coef_stderrs,
                            ci_lower = coef_ests - (qnorm(0.975) * coef_stderrs),
                            ci_upper = coef_ests + (qnorm(0.975) * coef_stderrs),
                            pvalue = 2 * pnorm(abs(coef_ests / coef_stderrs), lower.tail = FALSE))
results_table <- merge(results_table, target_snps, by.x = "id", by.y = "ID", all.x = TRUE, sort = FALSE)