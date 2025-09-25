#### fit issa's, generate personality specific spreads, and investigate covariates ####

rm(list = ls())

library(amt)
library(tidyverse)
library(lubridate)
library(spdep)
library(survival)
library(rasterVis)
library(corrplot)
library(pROC)   # For ROC analysis
library(dismo)  # For k-fold cross-validation (to split the data)
library(MASS)   # For logistic regression models (if used)
# Run git lfs pull to fetch large files before loading
system("git lfs pull")
# load data
load("data/implicit_data.rda")
# get movement data attributes
move_scaling<- implicit_used_control %>%
  group_by(pumaID) %>%
  summarize(mean_sl = mean(step_dist),
            sd_sl = sd(step_dist),
            mean_ta = mean(turning_angle),
            sd_ta = mean(turning_angle))
# and scale movement data per individual
to_model<- implicit_used_control %>%
  left_join(., move_scaling, by = "pumaID") %>%
  mutate(sl_scaled = (step_dist - mean_sl)/sd_sl,
         ta_scaled = (turning_angle - mean_ta)/sd_ta)

# check for correlation
# Visualize the correlation matrix
vars<- to_model %>% filter(case_ == 0) %>% dplyr::select(hd_150, 
                                                         prox_urbanedge,
                                                         slope, 
                                                         cover)
plot.new()
par(mfrow = c(1,1))
corrplot::corrplot(cor(vars, use = "complete.obs"), type = "upper", order = "hclust",
                   tl.col = "black", tl.srt = 45, method = "number") #plot correlation

# fit pop model to individuals
fits<- list()
coef_db<- data.frame(pumaID = unique(to_model$pumaID),
                     hd_150  = rep(NA, length(unique(to_model$pumaID))),
                     prox_urbanedge = rep(NA, length(unique(to_model$pumaID))),
                     slope = rep(NA, length(unique(to_model$pumaID))),
                     cover = rep(NA, length(unique(to_model$pumaID))),
                     landcover_simphuman = rep(NA, length(unique(to_model$pumaID))),
                     log_sl = rep(NA, length(unique(to_model$pumaID))),
                     step_dist = rep(NA, length(unique(to_model$pumaID))),
                     turning_angle = rep(NA, length(unique(to_model$pumaID))),
                     status = rep(NA, length(unique(to_model$pumaID)))
)

#### fit models ####
individuals<- split(to_model, to_model$pumaID)
for (i in 1:length(unique(to_model$pumaID))) {
  fit_df<- individuals[[i]]  
  fits[[i]]<- clogit(case_ ~ 
                       hd_150 + 
                       prox_urbanedge  +
                       slope + 
                       cover  + 
                       landcover_simp +
                       log_sl +
                       sl_scaled + 
                       ta_scaled + 
                       strata(end_gpsID),
                     fit_df, model = TRUE)
  coef_db[i, 2:9]<- unname(coef(fits[[i]]))
  coef_db$status[i]<- ifelse(individuals[[i]]$resident[1] == TRUE, "res", "disp")
}

# extract individual coefs
detailed_summary <- lapply(fits, function(model) {
  model_summary <- summary(model)
  data.frame(Coefficient = coef(model),
             StdError = sqrt(diag(vcov(model))))
}) %>% 
  as.data.frame() %>% 
  mutate(var = row.names(.)) %>% 
  pivot_longer(., names_to = "term", values_to = "coef", cols = c(1:length(colnames(.)) - 1)) %>% 
  mutate(pumaID = rep(unique(to_model$pumaID), each = 2, times = 8),
         term = gsub("\\.[0-9]+$", "", term)) %>%
  pivot_wider(., id_cols = c("var", "pumaID"), names_from = "term", 
              values_from = "coef") %>%
  mutate(upper = Coefficient + 1.96*StdError,
         lower = Coefficient - 1.96*StdError) %>%
  filter(var %in% c("hd_150", "prox_urbanedge", 
                    "slope", "cover",
                    "landcover_simphuman"
  )) %>%
  mutate(var = case_when(var == "cover" ~ "tree cover",
                         var == "hd_150" ~ "building density",
                         var == "slope" ~ "slope",
                         var == "prox_urbanedge" ~ "distance to urban center",
                         var == "landcover_simphuman" ~ "urban or ag landcover"
                         ))

#### cross validation procedure #####
ranks<- list()
kfoldSSF<- function(mod, k, nrepet, jitter, reproducible){
  ## Try to retrieve the data
  dt <- try(model.frame(mod), silent = TRUE)
  ## If it failed, stop and give a solution
  if (class(dt) == "try-error")
    stop("'model.frame' was unable to retrieve the data.",
         "Use 'model = TRUE' in the 'coxph' or 'clogit' call.")
  ## The first column is named 'srv' instead of 'Surv(faketime,
  ## case)'
  names(dt)[1] <- "srv"
  ## Which column is the strata?
  nstr <- attr(terms(mod), "specials")$strata
  ## Ugly regexp to extract and apply the strata variable name
  names(dt)[nstr] <- namestr <- sub("strata\\((.*)\\)", "\\1",
                                    names(dt)[nstr])
  ## If there is a cluster term...
  if (!is.null(attr(terms(mod), "specials")$cluster)) {
    ## Which column is the cluster?
    nclu <- attr(terms(mod), "specials")$cluster
    ## Ugly regexp to extract and apply the cluster variable name
    names(dt)[nclu] <- sub("cluster\\((.*)\\)", "\\1",
                           names(dt)[nclu])
  }
  ## Is it really a problem?
  ## ncase <- table(tapply(dt$srv[, 2], dt[, nstr], function(x) sum(x == 1)))
  ## if (any(names(ncase) == "0"))
  ##     stop(paste("Some stratas had no case.",
  ##       "It is likely that NAs were present in the variables for some cases."))
  ## Prepare the 'kfold', 'rd' and 'warn' objects
  kfold <- rd <- warn <- numeric(length = 100)
  
  ## The core of the kfold, each repetition
  for (i in 1:nrepet) {
    ## Create a 'set' column, which defaults to "train"
    dt$sets <- "train"
    ## Allows for reproducibility
    if (reproducible)
      set.seed(i)
    ## Sample the "test" data set
    dt$sets[dt[, namestr] %in% sample(unique(dt[, namestr]),
                                      length(unique(dt[, namestr]))/k)] <- "test"
    ## Update the regression using the training data
    reg <- update(mod, srv ~ ., data = subset(dt, sets ==
                                                "train"), model = TRUE)
    ## Extract the "test" data set
    dtest <- droplevels(subset(dt, sets == "test"))
    ## And compute the predictions associated to this data set
    ## using the training regression
    dtest$predall <- exp(predict(reg, type = "lp", newdata = dtest,
                                 reference = "sample"))
    ## In case of equality among predictions (e.g. categorical
    ## variable), add some noise to the predictions
    if (jitter) {
      ## Allows for reproducibility
      if (reproducible)
        set.seed(i)
      dtest$predall <- jitter(dtest$predall)
    }
    ## The function to compute the rank within a strata
    samplepred <- function(df) {
      ## Number of controls
      nrand <- sum(df$srv[, 2] == 0)
      ## Rank of the case (among case + controls)
      obs <- rank(df$predall)[df$srv[, 2] == 1]
      ## Rank of a random control (among controls only!)
      if (reproducible)
        set.seed(i)
      rand <- sample(rank(df$predall[df$srv[, 2] == 0]),
                     1)
      return(data.frame(obs = obs, rand = rand, nrand = nrand))
    }
    ## Compute the ranks for each strata and bind them together
    ranks <- do.call(rbind, by(dtest, dtest[, namestr], function(df) {
      # Only process strata with at least one case
      if (sum(df$srv[, 2] == 1) == 0) {
        return(NULL)  # Skip strata with no cases
      } else {
        samplepred(df)  # Process normally otherwise
      }
    }))
    ## Is there the same number of controls per strata?
    nrand <- unique(ranks$nrand)
    ## If no, use the greatest number of controls (and keep track
    ## of it)
    if (length(nrand) != 1) {
      nrand <- max(nrand)
      warn[i] <- 1
    }
    ## Compute the Spearman correlation on the ranks for the cases
    kfold[i] <- cor(1:(nrand+1), table(factor(ranks$obs,
                                              levels = 1:(nrand+1))), method = "spearman")
    ## Same for the random controls
    rd[i] <- cor(1:(nrand), table(factor(ranks$rand, levels = 1:(nrand))),
                 method = "spearman")
  }
  ## Create a data frame with the correlations and the type of value
  res <- data.frame(cor = c(kfold, rd), type = rep(c("obs",
                                                     "rand"), each = nrepet))
  ## Return the result data frame
  return(res)
}
for (cross in 1:length(fits)) {
  ranks[[cross]]<- kfoldSSF(fits[[cross]], k = 5, nrepet = 100, jitter = TRUE, reproducible = TRUE)
}
# Summarize performance for each model
summary_ranks <- lapply(ranks, function(res) {
  obs_cor <- res[res$type == "obs", "cor"]  # Extract observed correlations
  data.frame(
    mean_cor = mean(obs_cor),  # Mean correlation
    sd_cor = sd(obs_cor),      # Standard deviation
    min_cor = min(obs_cor),    # Minimum correlation
    max_cor = max(obs_cor)     # Maximum correlation
  )
})
# Combine results into a single data frame for easier comparison
summary_df <- do.call(rbind, summary_ranks)
rownames(summary_df) <- paste0("Individual_", seq_along(ranks))
summary_df
# Combine all results into one data frame for plotting
library(ggplot2)

ranks_combined <- do.call(rbind, lapply(seq_along(ranks), function(i) {
  data <- ranks[[i]]
  data$individual <- paste0("Individual_", i)
  data
}))

# Boxplot of observed correlations per individual
ggplot(ranks_combined[ranks_combined$type == "obs", ], aes(x = individual, y = cor)) +
  geom_boxplot() +
  labs(title = "Model Performance Per Individual",
       x = "Individual",
       y = "Spearman Correlation (Observed)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine observed correlations across all individuals
all_obs <- do.call(rbind, lapply(ranks, function(res) res[res$type == "obs", "cor"]))

# Calculate overall performance metrics
overall_summary <- data.frame(
  mean_obs = mean(all_obs),          # Mean observed correlation
  sd_obs = sd(all_obs),              # Standard deviation
  min_obs = min(all_obs),            # Minimum observed correlation
  max_obs = max(all_obs)             # Maximum observed correlation
)
overall_summary

#### calculate weighted average coefs ####
# Function to extract coefficients and standard errors from a model
extract_info <- function(model) {
  # Extract coefficients
  coefs <- coef(model)
  
  # Extract standard errors
  se <- sqrt(diag(vcov(model)))
  
  # Return a list of coefficients and standard errors
  return(list(coefs = coefs, se = se))
}
# Add sex information to the weights dataframe
sex_info <- to_model %>%
  mutate(sex = tolower(gsub('[0-9]+', "", pumaID))) %>%
  distinct(pumaID, sex) %>%
  dplyr::select(sex) %>%
  mutate(num = 1:n())
males<- sex_info[sex_info$sex == "m", "num"]
females<- sex_info[sex_info$sex == "f", "num"]
# Subset models based on indices
male_models<- list()
female_models<- list()
for (i in 1:length(males)) {
  male_models[[i]]<- fits[[males[i]]]
}
for (i in 1:length(females)) {
  female_models[[i]]<- fits[[females[i]]]
}
# Extract information from each model
model_info_m <- lapply(male_models, extract_info)
model_info_f<- lapply(female_models, extract_info)
# Get covariate names (assuming they are the same across all models)
covariate_names <- names(model_info_m[[1]]$coefs)

# Calculate weights: you can use 1/(se^2) or another method
# Here, we use inverse of variance (1/se^2): this means that pumas with more precision (smaller standard errors) are weighted heavier
weights_m_unnorm <- sapply(model_info_m, function(info) 1 / (info$se^2))
weights_f_unnorm <- sapply(model_info_f, function(info) 1 / (info$se^2))

# Normalize weights for each sex
weights_m_unnorm[5,23]<- 0 # dealing with one problematic animal; assigning a weight of 0
weights_m_unnorm[5,35]<- 0
weights_m<- weights_m_unnorm / rowSums(weights_m_unnorm, na.rm = TRUE) 
weights_f<- weights_f_unnorm / rowSums(weights_f_unnorm, na.rm = TRUE)
# Function to compute weighted quantiles
weighted_quantile <- function(values, weights, probs) {
  # Ensure that probs are within the range of 0 to 1
  if (any(probs < 0 | probs > 1)) {
    stop("probs should be between 0 and 1.")
  }
  
  # Order values and weights by values
  sorted_indices <- order(values)
  sorted_values <- values[sorted_indices]
  sorted_weights <- weights[sorted_indices]
  
  # Calculate cumulative weights and total weight
  cumulative_weights <- cumsum(sorted_weights)
  total_weight <- sum(sorted_weights)
  
  # Approximate quantiles based on cumulative weights
  quantiles <- approx(cumulative_weights / total_weight, sorted_values, xout = probs, rule = 2)$y
  
  return(quantiles)
}
# Compute weighted quantiles with custom probabilities
custom_probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# Calculate weighted averages, standard errors, and quantiles for each covariate
weighted_summary_m <- lapply(1:length(covariate_names), function(i) {
  coefs <- sapply(model_info_m, function(info) info$coefs[i])
  coefs[is.na(coefs)]<- 10000 # in case of any NA's, let function run then deal with later
  ses <- sapply(model_info_m, function(info) info$se[i])
  ses[is.na(ses)]<- 10000 # in case of any NA's, let function run then deal with later
  
  w <- weights_m[i, ]
  w[is.na(w)]<- 0
  
  weighted_avg <- sum(coefs * w)
  weighted_var <- sum((w * ses)^2)
  weighted_se <- sqrt(weighted_var)
  
  quantiles <- weighted_quantile(coefs, w, custom_probs)
  # Create named list for quantiles
  quantiles_named <- setNames(quantiles, paste0("p", custom_probs * 100))
  
  return(data.frame(
    cov = covariate_names[i],
    avg = weighted_avg,
    se = weighted_se,
    quantiles = quantiles_named
  ) %>%
    mutate(percentile = row.names(.)) %>%
    pivot_wider(., id_cols = c("cov", "avg", "se"), names_from = "percentile", values_from = "quantiles"))
})
weighted_summary_f <- lapply(1:length(covariate_names), function(i) {
  coefs <- sapply(model_info_f, function(info) info$coefs[i])
  ses <- sapply(model_info_f, function(info) info$se[i])
  w <- weights_f[i, ]
  
  weighted_avg <- sum(coefs * w)
  weighted_var <- sum((w * ses)^2)
  weighted_se <- sqrt(weighted_var)
  
  quantiles <- weighted_quantile(coefs, w, custom_probs)
  # Create named list for quantiles
  quantiles_named <- setNames(quantiles, paste0("p", custom_probs * 100))
  
  return(data.frame(
    cov = covariate_names[i],
    avg = weighted_avg,
    se = weighted_se,
    quantiles = quantiles_named
  ) %>%
    mutate(percentile = row.names(.)) %>%
    pivot_wider(., id_cols = c("cov", "avg", "se"), names_from = "percentile", values_from = "quantiles"))
})
# Convert lists to data frames
weighted_df_m <- do.call(rbind, weighted_summary_m)
weighted_df_f <- do.call(rbind, weighted_summary_f)
# Add sex column for clarity
weighted_df_m <- cbind(sex = "m", weighted_df_m)
weighted_df_f <- cbind(sex = "f", weighted_df_f)
# Combine the data frames
weighted_df <- rbind(weighted_df_m, weighted_df_f)
# Convert the list to a data frame for better presentation
weighted_df <- data.frame(var = rep(covariate_names, n = 2), weighted_df) %>%
  mutate(Covariate = case_when(var == "cover" ~ "tree cover",
                               var == "hd_150" ~ "building density",
                               var == "slope" ~ "slope",
                               var == "prox_urbanedge" ~ "distance to urban center",
                               var == "landcover_simphuman" ~ "urban or ag landcover",
                               var == "sl_scaled" ~ "step shape modifier",
                               var == "ta_scaled" ~ "direction modifier",
                               var == "log_sl" ~ "step rate modifier"))
# add max/min values
max_mins<- coef_db %>%
  pivot_longer(cols = hd_150:turning_angle,
               names_to = "variable", values_to = "value") %>%
  mutate(variable = gsub("\\.", ":", variable),
         sex = tolower(gsub("[^a-zA-Z]", "", pumaID))) %>% # changing out : for . for the two interaction terms
  group_by(variable, sex) %>%
  summarize(
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE)
  ) %>%
  ungroup()
weighted_df<- weighted_df %>% 
  left_join(., max_mins, by = c("var" = "variable", "sex"))
# arrange columns in order
# Desired column order
desired_order <- c("var", "sex", "min", "p10", 
                   "p20", "p30", "p40", "p50", "avg", 
                   "p60", "p70", "p80", "p90", "max",
                   "se", "Covariate")
# Rearrange the columns in the specified order
weighted_df <- weighted_df[, desired_order] %>%
  mutate(control = 0)
# Print the results
print(weighted_df)

#### save weighted data ####
save(coef_db, weighted_df, weights_m_unnorm, weights_f_unnorm, 
     file = "tmp/ind_coefs_and_weights.rda")

