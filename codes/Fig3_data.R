# This script specifies clinical data structure and contains functions used for Cox regression modeling to evaluate the association between clinical/microbial factors and the risk of developing a catheter-associated urinary tract infection (CAUTI).

library(survival)
library(coxme)

# Infection data was kindly provided by Dr. Chelsie Armbruster
# How Often Do Clinically Diagnosed Catheter-Associated Urinary Tract Infections in Nursing Homes Meet Standardized Criteria? Journal of the American Geriatrics Society 2016.

# Similar to the workflow for Fig 2, here we organize patient data individually, right censor when a CAUTI occurs, and run univariate and multivariate analyses to idenitfy polymicrobial interactions most strongly associated with a CAUTI with a specific organism (based on urine culture and clinical symptoms, but we do not have the resistance phenotype of the organism found in urine).

# Only include organisms that are found in 20 patients, and while CAUTI often is recurring, we right censor at the first event.

# Example dataset
cauti_data = readRDS("2020-07-04_org_CAUTI_mat.RData") 

names(cauti_data) # urine positive with each organism:
# Sa = Staphylococcus aureus; 
# Pa = Pseudomonas aeruginosa;
# Ec = Escherichia coli
# Enterococcus = Enterococcus
# Pm = Proteus mirabilis

colnames(cauti_data$Sa)
# start = enrolment date (set to 0)
# stop = one day after start
# Age = patient age
# Comorbidity.Score = Charlson comorbidity score (low = good, high = bad)
# sex = patient sex (0 = female; 1 = male)
# Dementia = patient dementia status (0 = no; 1 = yes)
# new_vpsms = physical self-maintenance score (0 = no assistance needed; 30 = completely dependent)
# urinepos = positive urine result with the organism of interest


# HR OF mono-colonization vs co-colonizatoin in the risk of CAUTI

# Only keeping the 4 most prevalent organisms retrieved from urine as outcomes (dropped P. aeruginosa)
org_by_urine_matrix = matrix(NA, nrow = 4, ncol = 4, 
                             dimnames = list(c("MRSA", "Escherichia coli", "VRE", "Proteus mirabilis"), 
                                             c("MRSA", "Escherichia coli", "VRE", "Proteus mirabilis")))

org_by_urine_matrix_pvalue = matrix(NA, nrow = 4, ncol = 4, 
                                    dimnames = list(paste0(c("MRSA", "Escherichia coli", "VRE", "Proteus mirabilis"), "_pvalue"), paste0(c("MRSA", "Escherichia coli", "VRE", "Proteus mirabilis"), "_pvalue")))

options(warn = 0)
gut_single_org_HR_by_org_in_urine = list()
org_ind = c("VRE", "MRSA", "Acinetobacter_baumanii", "Escherichia_coli", "Proteus_mirabilis","Pseudomonas_aeruginosa")

for (r in names(cauti_data)){print(r)
        
        if (r == "Sa"){o = "MRSA"};
        if (r == "Pa"){next};
        if (r == "Ec"){o = "Escherichia_coli"};
        if (r == "Enterococcus"){o ="VRE"};
        if (r == "Pm"){o = "Proteus_mirabilis"}

        org_urine_mat = matrix(NA, nrow = 3, ncol = 2, dimnames = list(c(o, "Others", "Both"), c("HR", "Wald p-value")))
        
        first_temp_df = cauti_data[[r]]
        
        o_col = which(colnames(first_temp_df) %in%  o)
        other_col = which(colnames(first_temp_df) %in%  org_ind[!org_ind %in% o])
        start_col = which(colnames(first_temp_df) %in% "start")
        stop_col = which(colnames(first_temp_df) %in% "stop")
        urine_col = which(colnames(first_temp_df) %in% "urinepos")
        psms_col = which(colnames(first_temp_df) %in% "new_vpsms")
        como_col = which(colnames(first_temp_df) %in% "Comorbidity.Score")
        age_col = which(colnames(first_temp_df) %in% "Age")
        trt_col = which(colnames(first_temp_df) %in% "Intervention.Site")
        sex_col = which(colnames(first_temp_df) %in% "sex")
        
        temp_df =  first_temp_df[,c(start_col, stop_col, urine_col, o_col, other_col, psms_col, como_col, age_col,sex_col, trt_col)]
        
        other_column = apply(temp_df[,org_ind[!org_ind %in% o]], 1, FUN = function(x){if (sum(x) > 0){1} else{0}})
        
        fac_column = substr(rownames(temp_df), 1, 1) # this needs to be changed depending on dataset
        id_column = substr(rownames(temp_df), 1, 4)

        temp_df = cbind(temp_df, other_column)
        
        org_type_vec = apply(temp_df[,c(o, "other_column")], 1, FUN = function(x){
                if(x[1] == 0 & x[2] == 0){org_type = 0}
                if(x[1] == 1 & x[2] == 0){org_type = 1}
                if(x[1] == 0 & x[2] == 1){org_type = 2}
                if(x[1] == 1 & x[2] == 1){org_type = 3}
                org_type
        })
        temp_df = cbind(temp_df, org_type_vec)
        
        org_mat = matrix(0, nrow = nrow(temp_df), ncol = 4, dimnames = list(NULL, c("none", "org", "others", "both")))
        
        # Create dummy variables indicating colonization patterns 
        # none: no colonization
        # org: colonization with the organism of interest 
        # others: colonization with any other organism 
        # both: colonization with organism of interest and at least another organism
        
        for (i in colnames(org_mat)){print(i)
                if (i == "none"){
                        org_vec_ind = which(org_type_vec == 0)
                        org_mat[org_vec_ind,i] = 1 
                }
                
                if (i == "org"){
                        org_vec_ind = which(org_type_vec == 1)
                        org_mat[org_vec_ind,i] = 1
                }
                
                if (i == "others"){
                        org_vec_ind = which(org_type_vec == 2)
                        org_mat[org_vec_ind,i] = 1
                }
                
                if (i == "both"){
                        org_vec_ind = which(org_type_vec == 3)
                        org_mat[org_vec_ind,i] = 1
                }
        }
        temp_df = cbind(temp_df, org_mat)
        
        if (sum(is.na(temp_df)) > 0){break}
        
        if(sum(table(temp_df$urinepos, temp_df$org_type_vec)[2,] == 0) > 0){next}
        
        mono_coxme_model = coxme(Surv(start, stop, urinepos) ~ temp_df$org + temp_df$others + temp_df$both + temp_df$new_vpsms + temp_df$sex + (1 | fac_column/id_column), data = temp_df)
        
        x = mono_coxme_model
        beta <- x$coefficients
        nvar <- length(beta)
        nfrail <- nrow(x$var) - nvar
        
        se <- sqrt(diag(x$var)[nfrail + 1:nvar])
        tmp <- cbind(beta, exp(beta), se, round(beta/se, 4), 
                     signif(1 - pchisq((beta/se)^2, 1), 4))
        
        coxme_model_pvalue = round(tmp[1:3,5],4)
        coxme_model_HR = round(tmp[1:3,2],4)
        
        options(scipen=999)
        
        org_urine_mat[rownames(org_urine_mat) %in% o,]  = c(coxme_model_HR[1], coxme_model_pvalue[1])
        org_urine_mat[rownames(org_urine_mat) %in% "Others",]= c(coxme_model_HR[2], coxme_model_pvalue[2])
        org_urine_mat[rownames(org_urine_mat) %in% "Both",]= c(coxme_model_HR[3], coxme_model_pvalue[3])
        
        gut_single_org_HR_by_org_in_urine[[paste0(o, " in urine")]] = org_urine_mat
        
        o2 = gsub("_", " ", o)
        org_by_urine_matrix[o2,o2] = org_urine_mat[o,"HR"]
        org_by_urine_matrix_pvalue[paste0(o2, "_pvalue"),paste0(o2, "_pvalue")] = org_urine_mat[o, "Wald p-value"]
        mono_outcome_list_for_print[[o]] = tmp
        
}
        
#HR OF CO-COLONIZATION
gut_cocolonization_org_HR_by_org_in_urine = list()
outcome_orgs = c("MRSA", "Escherichia_coli", "VRE", "Proteus_mirabilis")
org_pairs = t(combn(outcome_orgs, 2))
for (r in names(cauti_data)){print(r)
        
        if (r == "Sa"){o = "MRSA"};
        if (r == "Pa"){next};
        if (r == "Ec"){o = "Escherichia_coli"};
        if (r == "Enterococcus"){o ="VRE"};
        if (r == "Pm"){o = "Proteus_mirabilis"}
 
gut_cocolonization_org_pairs_HR_by_org_in_urine = list()   
outcome_list_for_print_by_org = list()
co_summary_table_by_org = list()

for (s in 1:nrow(org_pairs)){print(s)
        org1 = org_pairs[s,1]        
        org2 = org_pairs[s,2]        
        org_urine_mat = matrix(NA, nrow = 3, ncol = 2, dimnames = list(c(org1, org2, "Both"), c("HR", "Wald p-value")))
        
        first_temp_df = cauti_data[[r]]
        
        org1_col = which(colnames(first_temp_df) %in%  org1)
        org2_col = which(colnames(first_temp_df) %in%  org2)
        start_col = which(colnames(first_temp_df) %in% "start")
        stop_col = which(colnames(first_temp_df) %in% "stop")
        urine_col = which(colnames(first_temp_df) %in% "urinepos")
        psms_col = which(colnames(first_temp_df) %in% "new_vpsms")
        como_col = which(colnames(first_temp_df) %in% "Comorbidity.Score")
        age_col = which(colnames(first_temp_df) %in% "Age")
        trt_col = which(colnames(first_temp_df) %in% "Intervention.Site")
        sex_col = which(colnames(first_temp_df) %in% "sex")
        
        temp_df =  first_temp_df[,c(start_col, stop_col, urine_col, org1_col, org2_col, psms_col, como_col, age_col, sex_col, trt_col)]
        fac_column = substr(rownames(temp_df), 1, 1)
        id_column = substr(rownames(temp_df), 1, 4)

        org_type_vec = apply(temp_df[,c(org1, org2)], 1, FUN = function(x){
                if(x[1] == 0 & x[2] == 0){org_type = 0}
                if(x[1] == 1 & x[2] == 0){org_type = 1}
                if(x[1] == 0 & x[2] == 1){org_type = 2}
                if(x[1] == 1 & x[2] == 1){org_type = 3}
                org_type
        })
        temp_df = cbind(temp_df, org_type_vec, fac_column, id_column)
        
        org_mat = matrix(0, nrow = nrow(temp_df), ncol = 4, dimnames = list(NULL, c("none", "org1", "org2", "Both")))
        
        for (i in colnames(org_mat)){print(i)
                if (i == "none"){
                        org_vec_ind = which(org_type_vec == 0)
                        org_mat[org_vec_ind,i] = 1 
                }
                
                if (i == "org1"){
                        org_vec_ind = which(org_type_vec == 1)
                        org_mat[org_vec_ind,i] = 1
                }
                
                if (i == "org2"){
                        org_vec_ind = which(org_type_vec == 2)
                        org_mat[org_vec_ind,i] = 1
                }
                
                if (i == "Both"){
                        org_vec_ind = which(org_type_vec == 3)
                        org_mat[org_vec_ind,i] = 1
                }
        }
        temp_df = cbind(temp_df, org_mat)
        if(sum(is.na(temp_df)) > 0){break}
        
        if(sum(table(temp_df$urinepos, temp_df$org_type_vec)[2,] == 0) > 0){next}
        
        co_summary_table_by_org[[paste(org1, org2, sep = "/")]] = table(temp_df[temp_df$urinepos == 1, "org_type_vec"])
                
        co_coxme_model = coxme(Surv(start, stop, urinepos) ~temp_df$org1 + temp_df$org2 + temp_df$Both +      
                                       temp_df$new_vpsms +  temp_df$sex + (1 | fac_column/id_column),  data = temp_df)
        
        #EXTRACT VALUES FROM COXME_MODEL
        x = co_coxme_model
        beta <- x$coefficients
        nvar <- length(beta)
        nfrail <- nrow(x$var) - nvar
        
        se <- sqrt(diag(x$var)[nfrail + 1:nvar])
        tmp <- cbind(beta, exp(beta), se, round(beta/se, 4), 
                     signif(1 - pchisq((beta/se)^2, 1), 4))
        
        coxme_model_pvalue = round(tmp[1:3,5],4)
        coxme_model_HR = round(tmp[1:3,2],4)
        
        options(scipen=999)
        
        org_urine_mat[rownames(org_urine_mat) %in% org1,]  = c(coxme_model_HR[1], coxme_model_pvalue[1])
        org_urine_mat[rownames(org_urine_mat) %in% org2,]= c(coxme_model_HR[2], coxme_model_pvalue[2])
        org_urine_mat[rownames(org_urine_mat) %in% "Both",]= c(coxme_model_HR[3], coxme_model_pvalue[3])
        
        
        second_org = c(org1, org2)[!c(org1, org2) %in% o]
        
        if (length(second_org) == 2)
                
        {next}
        
        second_org2 = gsub("_", " ", second_org)
        o2 = gsub("_", " ", o)
        
        org_by_urine_matrix[second_org2,o2] = org_urine_mat["Both","HR"]
        org_by_urine_matrix_pvalue[paste0(second_org2, "_pvalue"),paste0(o2, "_pvalue")] = org_urine_mat["Both", "Wald p-value"]
        
        
        gut_cocolonization_org_pairs_HR_by_org_in_urine[[paste(org_pairs[s,1], org_pairs[s,2], sep = "/")]] = org_urine_mat
}
       
gut_cocolonization_org_HR_by_org_in_urine[[paste0(o, " in urine")]] = gut_cocolonization_org_pairs_HR_by_org_in_urine
}
        
# gut_cocolonization_org_HR_by_org_in_urine can be accessed by reading:
org_by_urine_matrix = readRDS("2020-07-04_org_cauti_coxph_HR.RData")
org_by_urine_matrix_pvalue = readRDS("2020-07-04_org_cauti_coxph_pval.RData")
org_by_urine_list = readRDS("2020-07-04_org_cauti_coxph_HR_list.RData") # same information as the matrices but in list

# Transform results into a matrix
outcome_network_matrix = matrix(NA, nrow = 8, ncol = 8 , 
                                dimnames = list(c("VRE colonization", "MRSA colonization", "P. mirabilis colonization", "E. coli colonization", "VRE CAUTI", "MRSA CAUTI", "P. mirabilis CAUTI", "E. coli CAUTI"), 
c("VRE colonization", "MRSA colonization", "P. mirabilis colonization", "E. coli colonization", 
 "VRE CAUTI", "MRSA CAUTI", "P. mirabilis CAUTI", "E. coli CAUTI")))

# In addition to associations between colonizatoin and CAUTI, we also want to overlay the interactions between microbial colonizations.

# Load Fig2 final matrix:
final_mat = readRDS("2020-07-04_org_abx_mat.RData")

# First fill interactions between microbial colonizations
for (i in 1:4){print(i)
        org = substr(rownames(outcome_network_matrix)[i], 1, nchar(rownames(outcome_network_matrix)[i]) - sum(nchar("colonization") + 1))
        if (org == "P. mirabilis"){org = "Proteus_mirabilis"}
        if (org == "E. coli"){org = "Escherichia_coli"}
        
        rfHR = final_mat[org,][which(final_mat[org,] > 0)]
        rf = names(which(final_mat[org,] > 0))
        rf_new = gsub("_", " ", rf)
        rf_new = sapply(rf_new, FUN = function(x){gsub(strsplit(x, " ")[[1]][1], paste0(substr(x,1,1), "."),x)})
        rf_new_colonize = paste0(rf_new, " colonization")
        rf_new_keep = which(rf_new_colonize %in% colnames(outcome_network_matrix))
        outcome_network_matrix[i, rf_new_colonize[rf_new_keep]] = rfHR[rf_new_keep]
        
}

# Next fill significant colonization-CAUTI associations
col_cauti_data_ind = which(org_by_urine_matrix_pvalue < 0.1, arr.ind = TRUE)

for (j in 1:nrow(col_cauti_data_ind)){print(j)
        outcome_org = colnames(org_by_urine_matrix)[col_cauti_data_ind[j,2]]
        rf_org = colnames(org_by_urine_matrix)[col_cauti_data_ind[j,1]]
        
        urine_spp_ind = c("MRSA", "Escherichia coli", "VRE", "Proteus mirabilis")
        
        outcome_ind = which(urine_spp_ind %in% outcome_org)
        
        cocol_ind = c(gsub(" ", "_", outcome_org), gsub(" ", "_", rf_org))
        cocol_ind = c(paste(cocol_ind, collapse  = "/"), paste(rev(cocol_ind), collapse = "/"))
        
        urine_data = org_by_urine_list[[outcome_ind]]
        urine_data = urine_data[names(urine_data) %in% cocol_ind]
        urine_data = urine_data[[1]]
        
        if (outcome_org == "Proteus mirabilis"){outcome_org = "P. mirabilis"}
        if (outcome_org == "Escherichia coli"){outcome_org = "E. coli"}
        
        rownames(urine_data)[rownames(urine_data) %in% "Escherichia_coli"] = "E. coli"
        rownames(urine_data)[rownames(urine_data) %in% "Proteus_mirabilis"] = "P. mirabilis"
        rownames(urine_data) = paste0(rownames(urine_data), " colonization")
        
        outcome_ind = grep(outcome_org, rownames(outcome_network_matrix)[5:8])
        outcome_network_matrix[5:8,c(rownames(urine_data)[1:2])][outcome_ind,] = urine_data[1:2, "HR"]
        
        
}

outcome_network_matrix[is.na(outcome_network_matrix)] = 0
colnames(outcome_network_matrix)[colnames(outcome_network_matrix) %in% "VRE CAUTI"] = "Enterococcus CAUTI"
rownames(outcome_network_matrix)[rownames(outcome_network_matrix) %in% "VRE CAUTI"] = "Enterococcus CAUTI"

# outcome_network_matrix is the final matrix containing HR for plotting Fig. 3 (saved as 2020-07-04_org_cauti_HR_mat.RData)
