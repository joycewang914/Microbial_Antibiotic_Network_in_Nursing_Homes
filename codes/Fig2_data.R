# This script specifies clinical data structure and contains functions used for Cox regression modeling to evaluate the association between clinical/microbial factors and colonization with different antibiotic-resistant organisms (ARO).

library(survival)
library(coxme)

# As a first pass, we loop through each combination of ARO colonization (response variable) and exposure to specific antibiotic or pre-existing ARO colonization (predictor variable)

# Example with ARO & antibiotics
# Prevalent ARO species
org_ind = c("VRE", "MRSA", "Acinetobacter_baumanii", "Escherichia_coli", "Proteus_mirabilis", "Pseudomonas_aeruginosa")

# Commonly used antibiotics
abx_ind = c("Aminoglycosides", "Carbapenems", "Glycopeptides", "Nitrofurans", "Nitroimidazoles", "Penicillins")

# Create all possible combinations of ARO/antibiotics and ARO/ARO
org_reps = rep(org_ind, each = length(abx_ind))
abx_reps = rep(abx_ind, times = length(org_reps)) 
org_abx_pairs = paste(org_reps, abx_reps, sep = "/")

# Sample clinical data
# Each patient is analyzed individually and data is arranged in a longitudinal format (follow-up days x clinical attributes). Each row indicates a day; once outcome (e.g. ARO colonization) occurs then data is right censored. e.g. patient 1 acquired VRE on day 121; patinet 2 never acquired an VRE during study so all follow-up days were included; patient 3 did not acquire an VRE during study.

patient1 = read.csv("2020-07-04_TIP_patient1_VRE_aminoglycoside.csv", header = T)
patient2 = read.csv("2020-07-04_TIP_patient2_VRE_aminoglycoside.csv", header = T)
patient3 = read.csv("2020-07-04_TIP_patient3_VRE_aminoglycoside.csv", header = T)
patient_list = list(pt1 = patient1, pt2 = patient2, pt3 = patient3)

# start = time of study enrolment (first row set to 0)
# stop = one day after start. 
# VRE = VRE culture result (0 = negative; 1 = positive)
# Aminoglycosides = Exposure to aminoglycosides within past 30 days (0 = no, 1 = yes)
# Comorbidity score = Charlson comorbidity score; usually constant over visits (low = good, high = bad)
# Age = patient age
# Sex = patient sex (0 = female; 1 = male)
# Dementia = patient dementia status (0 = no; 1 = yes)

# Processed individual patient data can be now combined in one large dataframe: 
total_pt_df = do.call(rbind, patient_list)

# Univariate analysis (use loop function to go through all combinations)
# Not adjusting for covariates or repeated measures yet
org = strsplit(org_abx_pairs[1], "/")[[1]][1]
abx = strsplit(org_abx_pairs[1], "/")[[1]][2]

abx_org_model = coxph(Surv(total_pt_df$start, total_pt_df$stop, total_pt_df[,org]) ~ total_pt_df[,abx])
summary_abx_org_model = summary(abx_org_model)
abx_org_model_coef = coef(summary_abx_org_model)
abx_org_hr = exp(abx_org_model_coef[1]) # Hazard ratio
abx_org_pval = abx_org_model_coef[4] # p-value

# Store coefficient and p-value in a matrix for each ARO/antibiotic combination as a matrix

# Load data frame containing p-value
# ARO colonization and pre-existing ARO colonization (Table S3)
second_org_as_outcome_pvalue = readRDS("2020-07-04_org1_org2_coxph_pval.RData")

# ARO colonization and antibiotics (Table S4)
abx_org_as_outcome_pvalue = readRDS("2020-07-04_abx_org_coxph_pval.RData")

# Extract antibiotics and microbial factors associated with ARO colonization in univariate analysis (p < 0.1)
abx_org_p0.1 = names(abx_org_as_outcome_pvalue[which(abx_org_as_outcome_pvalue < 0.1),])
second_org_p0.1 = names(second_org_as_outcome_pvalue[which(second_org_as_outcome_pvalue < 0.1), ])

# Unique microbial factors
unique_org_p0.1 = unique(c(unlist(lapply(strsplit(abx_org_p0.1, "/"), FUN = function(x){x[1]})), unlist(lapply(strsplit(second_org_p0.1, "/"), FUN = function(y){y[2]}))))

# Unique antibiotics
unique_abx_p0.1 = sub("_30", "", unique(unlist(lapply(strsplit(abx_org_p0.1, "/"), FUN = function(x){x[2]}))))

# For each outcome ARO, group potential risk factors in a list
secondary_org_variates_p0.1 = list()
for (u in unique_org_p0.1){print(u)
  
  temp_name = paste0(sub("_", ".", u),"_as_secondary_org_variates_p0.1")
  
  secondary_org_variates_p0.1[[temp_name]] = c(
    unlist(lapply(strsplit(abx_org_p0.1[grep(u, abx_org_p0.1)], "/"), FUN = function(x){x[2]})),
    unlist(lapply(strsplit(second_org_p0.1, "/"), FUN = function(y){if(y[2] == u){y[1]}}))) 
}

names(secondary_org_variates_p0.1) = unique_org_p0.1

# Multivariate analysis and control for repeated measures at patient and facilit levels.

# As with univariate analysis, now loop through each patient for each outcome ARO and censor when appropriate. Instead of selecting only one predictor at a time, include antibiotics/pre-existing ARO colonization that are identified to be associated in univariate analysis, as well as patient ID, facility ID, and additional factors of interest. 

# For example:
patient4 = readRDS("2020-07-04_TIP_patient4_VRE_covariates.RData")
patient5 = readRDS("2020-07-04_TIP_patient5_VRE_covariates.RData")
patient6 = readRDS("2020-07-04_TIP_patient6_VRE_covariates.RData")
patient_list = list(pt4 = patient4, pt5 = patient5, pt6 = patient6) # a list of patients with each ARO outcome in a sublist
patient_list = lapply(patient_list, FUN = function(x){lapply(x, FUN = function(y){apply(y, 1:2, as.numeric)})}) # Make sure data is in numeric form for modeling

# Perform Mixed Effects Cox Model
two_org_abx30_secondary_org_HR_p0.1 = list()
for (o in names(secondary_org_variates_p0.1)){print(o)
  
  pt_org_abx30_data_p0.1 = lapply(patient_list, FUN = function(x){x[[o]]})
  
  pt_org_abx30_data_start_stop_p0.1 = list() # add start and stop days for each outcome
  
  for (p in names(patient_list)){print(p)
    temp_data = pt_org_abx30_data_p0.1[[p]]
    
    if (is.null(temp_data)){next} # patient was skipped during processing because they were colonized with ARO of interest at enrollment thus no longer at risk.
    
    mat = matrix(NA, nrow  = nrow(temp_data), ncol = 2)
    mat[,1] = 0:(nrow(mat) - 1)
    mat[,2] = 1:nrow(mat)
    mat = cbind(mat, temp_data)
    rownames(mat) = paste0(p, ".", 1:nrow(mat))
    pt_org_abx30_data_start_stop_p0.1[[p]] = mat
  }
  
  all_pt_org_abx30_data_p0.1 = as.data.frame(do.call(rbind, pt_org_abx30_data_start_stop_p0.1))
  
  temp_df = all_pt_org_abx30_data_p0.1 
  colnames(temp_df) = c("start", "stop", colnames(temp_df)[!colnames(temp_df) %in% c("V1", "V2")]) 
  
  temp_df = as.data.frame(sapply(temp_df, FUN = function(x){as.numeric(as.character(x))}))
  rownames(temp_df) = rownames(all_pt_org_abx30_data_p0.1)
  fac_id = substr(rownames(temp_df), 1,1) # fac_id (facility ID) and pt_id (patient ID) need to be modified depending on clinical data structure.
  pt_id = substr(rownames(temp_df), 1,4) # 
  
  temp_var = colnames(temp_df)[!colnames(temp_df) %in% c("start", "stop",o)]
  
  temp_var_list = c()
  
  for (t in temp_var){
    print(t)
    temp_var_paste = paste("temp_df[ ,", t, "]", sep = "'")
    temp_var_list = c(temp_var_list, temp_var_paste)
  }
  
  temp_var_list = c(temp_var_list, "(1 | fac_id/pt_id)")
  
  coxme_formula = paste("Surv(start, stop, temp_df[,o]) ~", paste(temp_var_list, collapse = " + "))
  coxme_formula = as.formula(coxme_formula)
  
  #all variables
  
  coxme_model = coxme(coxme_formula, data = temp_df)
  
  x = coxme_model
  beta <- x$coefficients
  nvar <- length(beta)
  nfrail <- nrow(x$var) - nvar
  
  se <- sqrt(diag(x$var)[nfrail + 1:nvar])
  tmp <- cbind(beta, exp(beta), se, round(beta/se, 4), 
               signif(1 - pchisq((beta/se)^2, 1), 4))
  
  
  two_org_abx30_secondary_org_HR_p0.1[[o]] = tmp
  
}

# Last, create a symmetric matrix summarizing modeling results 
# Only include hazard ratios with a p < 0.05 in multivariate analysis.
org_abx30_graph_mat_p0.05 = matrix(data = 0, 
                                   nrow = sum(length(unique_org_p0.1), 
                                              length(unique_abx_p0.1)), 
                                   ncol = sum(length(unique_org_p0.1), 
                                              length(unique_abx_p0.1)),
                                   dimnames = list(c(c(unique_org_p0.1, unique_abx_p0.1)), 
                                                   c(unique_org_p0.1, unique_abx_p0.1)))

for (o in names(two_org_abx30_secondary_org_HR_p0.1)){print(o)
  
  coef_HR = two_org_abx30_secondary_org_HR_p0.1[[o]][grep("sex", rownames(two_org_abx30_secondary_org_HR_p0.1[[o]]), invert = T),]
  
  
  sig_p0.05 = which(coef_HR[,5] < 0.05)
  sig_covariates_p0.05 = rownames(coef_HR)[sig_p0.05]
  
  sig_covariates_p0.05 = unlist(lapply(strsplit(sig_covariates_p0.05, "\"", fixed = TRUE), FUN = function(x){x[2]}))
  sig_covariates_p0.05 = sub("_30", "", sig_covariates_p0.05)
  org_abx30_graph_mat_p0.05[o, sig_covariates_p0.05] = coef_HR[coef_HR[,5] < 0.05,2]
  
}

org_abx30_graph_mat_p0.05_row_keep = names(which(apply(org_abx30_graph_mat_p0.05, 1, sum) > 0))
org_abx30_graph_mat_p0.05_col_keep = names(which(apply(org_abx30_graph_mat_p0.05, 2, sum) > 0))
org_abx30_graph_mat_p0.05_dimnames_keep = unique(union(org_abx30_graph_mat_p0.05_row_keep, org_abx30_graph_mat_p0.05_col_keep)) #USE UNIQUE TO AVOID DUPLICATED ROWS/COLUMNS

org_abx30_graph_mat_p0.05 = org_abx30_graph_mat_p0.05[org_abx30_graph_mat_p0.05_dimnames_keep, org_abx30_graph_mat_p0.05_dimnames_keep]
org_abx30_graph_mat_p0.05 = round(org_abx30_graph_mat_p0.05, 1)

# org_abx_30_graph_mat_p0.05 is the final matrix containing HR for plotting Fig. 2 (saved as 2020-07-04_org_abx_mat.RData)