#' ---
#' title: "Longitudinal study of risk and protective factors for gaming disorder"
#' author: "Ivan Ropovik"
#' date: "`r Sys.Date()`"
#' output:
#'    html_document:
#'       toc: true
#'       toc_float: true
#'       code_folding: show
#'       fig_retina: 2
#' always_allow_html: yes
#' ---
#+ setup, include = FALSE
knitr::opts_chunk$set(echo=FALSE, warning = FALSE, fig.align = "center")
rm(list = ls())

# TODO Qs
# Skus spravit: Latent Class Analysis or Growth Mixture Modeling: Are there distinct subgroups of individuals with different trajectories of IGD symptoms over time? What characterizes these different subgroups? Try no covariates. 
# Zatial nevieme: Cross-Lagged Panel Analysis: How do the variables influence each other over time? For example, does impulsivity at wave 1 predict IGD symptoms at wave 2, controlling for IGD symptoms at wave 1?
# Alternative mode pri LGCM aj ked tryCatch produkuje warning, nie len error
# Check velke kovariancie pri b-ss networku

# self-assessment von z prediktorov ale ako operacionalizacia IGD
# gdt + igds9sf ako primarne, self-assessment ako supplement

# Define fitmeasures
fitMeasures <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr")

# Define the number of bootstrap samples
nBoot <- 200

# Source scripts
source("functions.R")
source("dataWrangling.R")

# Variable names
predictors <- c("adhd_Latent", "net_addiction_Latent", "soc_media_addiction_Latent", 
                "impulsivity_Latent", "stress_Latent", "aggression_hostility_Latent", 
                "escape_Latent", "social_support_Latent", "selfassessment_Latent", 
                "anxiety_Latent", "depression_Latent", "social_anxiety_Latent")

# Create a list of symptom variable names
symptoms <- paste0("igd9sf_", 1:9)

# CFA results -------------------------------------------------------------

# Extract the fit measures from the results list
fitMeasuresList <- lapply(resultsList, function(results) {
  lapply(results, `[[`, "fitMeasures")
})

# Remove NA elements from the fitMeasuresList
fitMeasuresList <- removeNAElements(fitMeasuresList)

# Unnest and combine all data frames into one
summaryTable <- do.call(rbind, lapply(names(fitMeasuresList), function(wave) {
  do.call(rbind, lapply(names(fitMeasuresList[[wave]]), function(variable) {
    data.frame(fitMeasure = names(fitMeasuresList[[wave]][[variable]]),
               value = unlist(fitMeasuresList[[wave]][[variable]]),
               variable = variable,
               wave = wave)
  }))
}))

# Reshape the table to wide format
summaryTableWide <- summaryTable %>%
  spread(key = fitMeasure, value = value)

# Print the wide format summary table
round_df(summaryTableWide, 3)

# Reliability -------------------------------------------------------------
reliabilities <- lapply(dataList, computeReliability, variableList = variableList[-c(11:14)])
reliabilities_df <- do.call(rbind, lapply(reliabilities, data.frame))

# Create a list to store the two dataframes
reliabilities_list <- list()

# Split the dataframe into two: one for alpha and one for omega
reliabilities_list$alpha <- t(reliabilities_df[, grep("cronbachAlpha", names(reliabilities_df))])
reliabilities_list$omega <- t(reliabilities_df[, grep("omega", names(reliabilities_df))])

print(round_df(as.data.frame(reliabilities_list), 3))

# Longitudinal Changes ----------------------------------------------------
describe(dat[, predictors])

###
fitResults <- lapply(predictors, function(var) {
  tryCatch({
    igdModel <- paste(' 
      # latent variable
      iIgd =~ 1*igd9sf_Latent + 1*igd9sf_Latent.w2 + 1*igd9sf_Latent.w3
      sIgd =~ 0*igd9sf_Latent + 1*igd9sf_Latent.w2 + 2*igd9sf_Latent.w3
      
      iPredictor =~ 1*',var,' + 1*',var,'.w2 + 1*',var,'.w3
      sPredictor =~ 0*',var,' + 1*',var,'.w2 + 2*',var,'.w3
      
      iIgd ~ iPredictor
      sIgd ~ iPredictor + sPredictor
      
      iIgd ~~ sIgd
      iPredictor ~~ sPredictor
      
      # lag1 autocorrelations
      # igd9sf_Latent ~~ igd9sf_Latent.w2
      # igd9sf_Latent.w2 ~~ igd9sf_Latent.w3
      # ',var,' ~~ ',var,'.w2
      # ',var,'.w2 ~~ ',var,'.w3
    
    ', collapse = "")
    
    fitIGD <- growth(igdModel, data = dat, missing = "fiml", mimic = "Mplus", bootstrap = nBoot, optim.method = "BFGS")
    
    list(
      summary = summary(fitIGD, standardized = TRUE),
      fitMeasures = fitmeasures(fitIGD, fitMeasures),
      parameterEst = parameterestimates(fitIGD, standardized = TRUE)
    )
  }, error = function(e) {
    message(paste("Error in model with predictor: ", var, ". Trying alternative model..."))
    
    tryCatch({
      # alternative model
      igdModel <- paste(' 
        # latent variable
        iIgd =~ 1*igd9sf_Latent + 1*igd9sf_Latent.w2 + 1*igd9sf_Latent.w3
        sIgd =~ 0*igd9sf_Latent + 1*igd9sf_Latent.w2 + 2*igd9sf_Latent.w3
        
        iPredictor =~ 1*',var,' + 1*',var,'.w2 + 1*',var,'.w3
        sPredictor =~ 0*',var,' + 1*',var,'.w2 + 2*',var,'.w3
        
        iIgd ~ iPredictor
        sIgd ~ iPredictor + sPredictor
        
        iIgd ~~ sIgd
        iPredictor ~~ sPredictor
      ', collapse = "")
      
      fitIGD <- growth(igdModel, data = dat, missing = "fiml", mimic = "Mplus", bootstrap = nBoot, 
                       optim.method = "nlminb")
      
      list(
        summary = summary(fitIGD, standardized = TRUE),
        fitMeasures = fitmeasures(fitIGD, fitMeasures),
        parameterEst = parameterestimates(fitIGD, standardized = TRUE)
      )
    }, error = function(e) {
      message(paste("Error in alternative model with predictor: ", var))
      message("Error message: ", e$message)
      return(NULL) # return NULL for this model if there's an error
    })
  })
})
names(fitResults) <- predictors

summaryTable <- map_dfr(fitResults, function(result) {
  if (is.null(result) || !is.list(result)) {
    return(setNames(data.frame(rep(NA, 13)), c("chiSq", "df", "pvalue", "CFI", "TLI", "RMSEA", "iIgd_iPred", "sIgd_iPred", "sIgd_sPred", "iIgd_iPredP", "sIgd_iPredP", "sIgd_sPredP", "iIgd_iPredP_adj")))
  }
  
  fitIndices <- result$fitMeasures
  parameters <- result$parameterEst
  
  iIgd_iPredP <- parameters$pvalue[parameters$lhs == "iIgd" & parameters$rhs == "iPredictor"]
  sIgd_iPredP <- parameters$pvalue[parameters$lhs == "sIgd" & parameters$rhs == "iPredictor"]
  sIgd_sPredP <- ifelse(length(parameters$pvalue[parameters$lhs == "sIgd" & parameters$rhs == "sPredictor"]) > 0, parameters$pvalue[parameters$lhs == "sIgd" & parameters$rhs == "sPredictor"], NA)
  
  data.frame(
    chiSq = fitIndices["chisq"],
    df = fitIndices["df"],
    pvalue = fitIndices["pvalue"],
    CFI = fitIndices["cfi"],
    TLI = fitIndices["tli"],
    RMSEA = fitIndices["rmsea"],
    iIgd_iPred = parameters$std.all[parameters$lhs == "iIgd" & parameters$rhs == "iPredictor"],
    sIgd_iPred = parameters$std.all[parameters$lhs == "sIgd" & parameters$rhs == "iPredictor"],
    sIgd_sPred = ifelse(length(parameters$std.all[parameters$lhs == "sIgd" & parameters$rhs == "sPredictor"]) > 0, parameters$std.all[parameters$lhs == "sIgd" & parameters$rhs == "sPredictor"], NA),
    iIgd_iPredP = iIgd_iPredP,
    sIgd_iPredP = sIgd_iPredP,
    sIgd_sPredP = sIgd_sPredP,
    iIgd_iPredP_adj = p.adjust(iIgd_iPredP, method = "holm", n = length(result)),    # new adjusted p-value
    sIgd_iPredP_adj = p.adjust(sIgd_iPredP, method = "holm", n = length(result)),    # new adjusted p-value
    sIgd_sPredP_adj = p.adjust(sIgd_sPredP, method = "holm", n = length(result))     # new adjusted p-value
  )
}, .id = "Variable") %>% remove_rownames()

# Add the variable names from the original list to the summary table
summaryTable$Variable <- names(fitResults)

# Results table
round_df(summaryTable, 3)

#############
t1PredResults <- lapply(predictors, function(var) {
  igdModel <- paste(' 
      # latent variable
      iIgd =~ 1*igd9sf_Latent + 1*igd9sf_Latent.w2 + 1*igd9sf_Latent.w3
      sIgd =~ 0*igd9sf_Latent + 1*igd9sf_Latent.w2 + 2*igd9sf_Latent.w3
      
      iIgd ~ ',var,'
      sIgd ~ ',var,'
      
      iIgd ~~ sIgd
      
      # lag1 autocorrelations
      # igd9sf_Latent ~~ igd9sf_Latent.w2
      # igd9sf_Latent.w2 ~~ igd9sf_Latent.w3
    ', collapse = "")
  
  fitIGD <- growth(igdModel, data = dat, missing = "fiml", mimic = "Mplus", bootstrap = nBoot, optim.method = "BFGS")
  
  list(
    summary = summary(fitIGD, standardized = TRUE),
    fitMeasures = fitmeasures(fitIGD, fitMeasures),
    parameterEst = parameterestimates(fitIGD, standardized = TRUE)
  )
})
names(t1PredResults) <- predictors

t1PredSummaryTable <- map_dfr(t1PredResults, function(result) {
  if (is.null(result) || !is.list(result)) {
    return(setNames(data.frame(rep(NA, 13)), c("chiSq", "df", "pvalue", "CFI", "TLI", "RMSEA", "iIgd_pred", "sIgd_pred", "iIgd_predP", "sIgd_predP", "iIgd_predP_adj", "sIgd_predP_adj")))
  }
  
  fitIndices <- result$fitMeasures
  parameters <- result$parameterEst
  
  iIgd_predP <- parameters$pvalue[parameters$lhs == "iIgd" & parameters$op == "~"]
  sIgd_predP <- parameters$pvalue[parameters$lhs == "sIgd" & parameters$op == "~"]
  
  data.frame(
    chiSq = fitIndices["chisq"],
    df = fitIndices["df"],
    pvalue = fitIndices["pvalue"],
    CFI = fitIndices["cfi"],
    TLI = fitIndices["tli"],
    RMSEA = fitIndices["rmsea"],
    iIgd_pred = parameters$std.all[parameters$lhs == "iIgd" & parameters$op == "~"],
    sIgd_pred = parameters$std.all[parameters$lhs == "sIgd" & parameters$op == "~"],
    iIgd_predP = iIgd_predP,
    sIgd_predP = sIgd_predP,
    iIgd_predP_adj = p.adjust(iIgd_predP, method = "holm", n = length(result)),   
    sIgd_predP_adj = p.adjust(sIgd_predP, method = "holm", n = length(result))
  )
}, .id = "Variable") %>% remove_rownames()

# Add the variable names from the original list to the summary table
t1PredSummaryTable$Variable <- names(t1PredResults)

# Results table
round_df(t1PredSummaryTable, 3)

# Development of symptoms over time ---------------------------------------

# Initialize an empty list to store the results
resultsSymptoms <- list()

# Fit separate mixed-effects models for each symptom
for (symptom in symptoms) {
  # Create formulas for the null and full models
  null_formula <- as.formula(paste0(symptom, " ~ 1 + (1 + wave|ID)"))
  full_formula <- as.formula(paste0(symptom, " ~ wave + (1 + wave|ID)"))
  
  # Fit the null model
  null_model <- lmer(null_formula, data = datLong, REML = T)
  
  # Fit the full model
  full_model <- lmer(full_formula, data = datLong, REML = T)
  
  # Conduct the LR test
  lr_test <- anova(null_model, full_model)
  
  # Extract the wave coefficient and its standard error from the full model
  wave_coef <- fixef(full_model)["wave"]
  wave_se <- sqrt(diag(vcov(full_model))["wave"])
  
  # Extract the random effects
  random_effects <- as.data.frame(VarCorr(full_model)$ID)
  
  # Extract the standard deviations of the random effects
  intercept_sd <- sqrt(random_effects[1, 1])
  wave_sd <- sqrt(random_effects[2, 2])
  
  #Extract the residual SD
  residual_sd <- attr(VarCorr(full_model), "sc")
  
  # Extract the correlation between the random intercepts and slopes
  intercept_wave_corr <- random_effects[2, 1] / (intercept_sd * wave_sd)
  
  # Create a data frame row with the results
  symptom_results <- data.frame(
    symptom = symptom,
    wave_coef = wave_coef,
    wave_se = wave_se,
    intercept_sd = intercept_sd,
    wave_sd = wave_sd,
    residual_sd = residual_sd,
    intercept_wave_corr = intercept_wave_corr,
    lr_test_p_value = lr_test$`Pr(>Chisq)`[2]
  )
  
  # Store the model, its summary, and the results in a list
  symptom_list <- list(
    null_model = null_model,
    full_model = full_model,
    null_summary = summary(null_model),
    full_summary = summary(full_model),
    symptom_results = symptom_results
  )
  
  # Add the list to the results
  resultsSymptoms[[symptom]] <- symptom_list
}

# Create and format a summary table with the results for all symptoms
summaryTableSymptoms <- do.call(rbind, lapply(resultsSymptoms, function(x) x$symptom_results)) %>% remove_rownames()

# Rename the columns and rows
names(summaryTableSymptoms) <- c("Symptom", "Wave Coefficient", "SE", "Intercept SD", "Wave SD", "Residual SD", "Correlation", "p-value")
summaryTableSymptoms$Symptom <- gsub("igd9sf_", "IGDS9-SF ", summaryTableSymptoms$Symptom)

# Print the table 
round_df(summaryTableSymptoms, 3)
# The negative correlations suggest that individuals who start with higher symptom levels tend to experience a smaller change or improvement in their symptoms over time.

# Wave Coefficient: the estimated fixed effect of wave on the symptom, indicating the average change in the symptom per unit increase in wave.
# Standard Error: the standard error of the wave coefficient, providing a measure of its uncertainty.
# Intercept SD and Wave SD: the standard deviations of the random intercepts and slopes, respectively, indicating the variability in these effects across participants.
# Intercept-Wave Correlation: the correlation between the random intercepts and slopes, indicating the extent to which participants with higher baseline symptom levels also tend to show larger wave effects.
# LR Test p-value: the p-value from the likelihood ratio test comparing the full model (including wave as a fixed effect and random slope) to the null model (excluding wave), providing evidence for whether wave has a significant effect on the symptom.

# Plot the fitted values
# Extract the fitted values from each model
fitted_values <- data.frame()
for (symptom in symptoms) {
  # Get the data used in the model
  model_data <- getME(resultsSymptoms[[symptom]][["full_model"]], "X")
  
  # Extract the fitted values
  symptom_fitted_values <- data.frame(
    wave = as.data.frame(model_data)$wave,
    symptom = rep(gsub("igd9sf_", "IGDS9-SF ", symptom), nrow(model_data)),
    fitted_value = fitted(resultsSymptoms[[symptom]][["full_model"]])
  )
  
  # Add the fitted values to the data frame
  fitted_values <- rbind(fitted_values, symptom_fitted_values)
}

# Calculate the mean fitted value for each wave and symptom
mean_fitted_values <- aggregate(fitted_value ~ wave + symptom, data = fitted_values, FUN = mean)

# Create the plot
ggplot(mean_fitted_values, aes(x = wave, y = fitted_value, color = symptom)) +
  geom_line() +
  scale_x_discrete(breaks = c("1", "2", "3"), expand = c(0,0)) + 
  scale_color_brewer(palette = "Set1") +
  labs(x = "Wave", y = "Mean Fitted Value", color = "Symptom") +
  theme_minimal() +
  ggtitle("Mean Fitted Symptom Value Over Time")

# Development of risk factors over time -----------------------------------
# Initialize an empty list to store the results
resultsPredictors <- list()

# Fit separate mixed-effects models for each risk factor
for (predictor in predictors) {
  # Create formulas for the null and full models
  null_formula <- as.formula(paste0(predictor, " ~ 1 + (1 + wave|ID)"))
  full_formula <- as.formula(paste0(predictor, " ~ wave + (1 + wave|ID)"))
  
  # Fit the null model
  null_model <- lmer(null_formula, data = datLong, REML = T)
  
  # Fit the full model
  full_model <- lmer(full_formula, data = datLong, REML = T)
  
  # Conduct the LR test
  lr_test <- anova(null_model, full_model)
  
  # Extract the wave coefficient and its standard error from the full model
  wave_coef <- fixef(full_model)["wave"]
  wave_se <- sqrt(diag(vcov(full_model))["wave"])
  
  # Extract the random effects
  random_effects <- as.data.frame(VarCorr(full_model)$ID)
  
  # Extract the standard deviations of the random effects
  intercept_sd <- sqrt(random_effects[1, 1])
  wave_sd <- sqrt(random_effects[2, 2])
  
  #Extract the residual SD
  residual_sd <- attr(VarCorr(full_model), "sc")
  
  # Extract the correlation between the random intercepts and slopes
  intercept_wave_corr <- random_effects[2, 1] / (intercept_sd * wave_sd)
  
  # Create a data frame row with the results
  predictor_results <- data.frame(
    predictor = predictor,
    wave_coef = wave_coef,
    wave_se = wave_se,
    intercept_sd = intercept_sd,
    wave_sd = wave_sd,
    residual_sd = residual_sd,
    intercept_wave_corr = intercept_wave_corr,
    lr_test_p_value = lr_test$`Pr(>Chisq)`[2]
  )
  
  # Store the model, its summary, and the results in a list
  predictor_list <- list(
    null_model = null_model,
    full_model = full_model,
    null_summary = summary(null_model),
    full_summary = summary(full_model),
    predictor_results = predictor_results
  )
  
  # Add the list to the results
  resultsPredictors[[predictor]] <- predictor_list
}

# Create and format a summary table with the results for all risk factors
summaryTablePredictors <- do.call(rbind, lapply(resultsPredictors, function(x) x$predictor_results)) %>% remove_rownames()

# Rename the columns and rows
names(summaryTablePredictors) <- c("Predictor", "Wave Coefficient", "SE", "Intercept SD", "Wave SD", "Residual SD", "Correlation", "p-value")
summaryTablePredictors$Predictor <- gsub("igd9sf_", "IGDS9-SF ", summaryTablePredictors$Predictor)

# Print the table 
round_df(summaryTablePredictors, 3)

# TODO Check code
# Plot the fitted values
# Extract the fitted values from each model
predFitted <- data.frame()
for (predictor in predictors) {
  # Get the data used in the model
  model_data <- getME(resultsPredictors[[predictor]][["full_model"]], "X")

  # Extract the fitted values
  predictorFittedValues <- data.frame(
    wave = as.data.frame(model_data)$wave,
    predictor = rep(gsub("igd9sf_", "IGDS9-SF ", predictor), nrow(model_data)),
    fitted_value = fitted(resultsPredictors[[predictor]][["full_model"]])
  )

  # Add the fitted values to the data frame
  predFitted <- rbind(predFitted, predictorFittedValues)
}

# Calculate the mean fitted value for each wave and symptom
predMeanFitted <- aggregate(fitted_value ~ wave + predictor, data = predFitted, FUN = mean)

# Create the plot
ggplot(predMeanFitted, aes(x = wave, y = fitted_value, color = predictor)) +
  geom_line() +
  scale_x_discrete(breaks = c("1", "2", "3"), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  #scale_color_brewer(palette = "Set1") +
  labs(x = "Wave", y = "Mean Fitted Value", color = "Risk factor") +
  theme_minimal() +
  ggtitle("Mean Fitted Risk Factor Value Over Time")

# Panel VAR network model -------------------------------------------------

# Form maximum likelihood covariance matrix
datVAR <- dat %>% select(starts_with("igd9sf_") & !contains("Latent"))
nVAR <- nrow(datVAR)
covVAR <- (datVAR - 1)/datVAR * cov(datVAR, use = "pairwise.complete.obs")
meansVAR <- colMeans(datVAR,  na.rm = TRUE)

# Define design
designVAR <- matrix(dat %>% select(starts_with("igd9sf_") & !contains("Latent")) %>% colnames(),nrow = 9, ncol = 3)

model <- panelgvar(data = datVAR, means = meansVAR, vars = designVAR,
                   identification = c("loadings","variance"), 
                   beta = "full",
                   covtype = "ML",
                   estimator = "FIML")

# Unprunned model
model <- model %>% runmodel(bounded = FALSE) %>% print

# Prunned model
# Pruning parameters
alpha <- 0.5
adjust <- "none"
model2 <- model %>% prune(alpha = alpha, adjust = adjust, recursive = FALSE) #%>% modelsearch(verbose = TRUE)
# saveRDS(resultsVAR, file ="resultsVAR.RDS")

# Compare:
(comp <- compare(
  original = model,
  pruned = model2
))

# Model fit
fitModel <- model %>% fit
fitModel2 <- model2 %>% fit
model %>% parameters
model2 %>% parameters

# Fix the problematic parameter to zero:
# model <- model %>% fixpar()
saveRDS(model, file="model.RDS")

# Get matrices
temporal <- temporalCov <- getmatrix(model, "PDC")
contemporaneous <- contemporaneousCov <- getmatrix(model, "omega_zeta_within")
between <- betweenCov <- getmatrix(model, "omega_zeta_between")

# Labels
labels <- paste0("igd9sf_", 1:9)

# Plot networks:
layout(t(1:3))
qgraph(temporal, layout = "circle", labels = labels,
       title = "Temporal", shape = "rectangle", 
       #vsize = 30, vsize2 = 20,
       #mar = rep(8,4), asize = 7,
       theme = "colorblind")
box("figure")

qgraph(contemporaneous, labels = labels,
       title = "Contemporaneous",shape = "rectangle", 
       #vsize = 30, vsize2 = 20,
       #mar = rep(8,4),
       theme = "colorblind")
box("figure")

qgraph(between, layout = "circle", labels = labels,
       title = "Between-persons",shape = "rectangle", 
       #vsize = 30, vsize2 = 20,
       #mar = rep(8,4),
       theme = "colorblind")
box("figure")

# Parameter estimates
rownames(temporalCov) <- colnames(temporalCov) <- rownames(contemporaneousCov) <- colnames(contemporaneousCov) <- rownames(betweenCov) <- colnames(betweenCov) <- labels
# Temporal
temporalCor <- cov2cor(getmatrix(model, "PDC"))
temporalCov[lower.tri(temporalCov)] <- temporalCor[lower.tri(temporalCor)]
round_df(as.data.frame(temporalCov), 3)
# Contemporaneous
contemporaneousCor <- cov2cor(getmatrix(model, "sigma_zeta_within"))
contemporaneousCov[lower.tri(contemporaneousCov)] <- contemporaneousCor[lower.tri(contemporaneousCor)]
round_df(as.data.frame(contemporaneousCov), 3)
# Between
betweenCor <- cov2cor(getmatrix(model, "omega_zeta_between"))
betweenCov[lower.tri(betweenCov)] <- betweenCor[lower.tri(betweenCor)]
round_df(as.data.frame(betweenCov), 3)


##########################################
# Test parts
# Select relevant columns
dat_subset <- dat %>% select(aggression_hostility_1, aggression_hostility_2, aggression_hostility_3, aggression_hostility_4, aggression_hostility_5, aggression_hostility_6)

# Define the CFA model
cfa_model <- '
  aggression_hostility =~ aggression_hostility_1 + aggression_hostility_2 + aggression_hostility_3
                      + aggression_hostility_4 + aggression_hostility_5 + aggression_hostility_6
'

# Fit the model
cfa_fit <- cfa(cfa_model, data = dat_subset)

# Print summary of the model fit
summary(cfa_fit)
fitmeasures(cfa_fit, fitMeasures)
modindices(cfa_fit)

#########################################
# Growth curve model for igd9sf_Latent
igdModel <- ' 
  # latent variable
  iIgd =~ 1*igd9sf_Latent + 1*igd9sf_Latent.w2 + 1*igd9sf_Latent.w3
  sIgd =~ 0*igd9sf_Latent + 1*igd9sf_Latent.w2 + 2*igd9sf_Latent.w3
  
  iPredictor =~ 1*stress_Latent + 1*stress_Latent.w2 + 1*stress_Latent.w3
  sPredictor =~ 0*stress_Latent + 1*stress_Latent.w2 + 2*stress_Latent.w3
  
  iIgd ~ iPredictor
  sIgd ~ iPredictor + sPredictor
  
  iIgd ~~ sIgd
  iPredictor ~~ sPredictor
'

fitIGD <- growth(igdModel, data = dat, missing = "fiml", mimic = "Mplus", bootstrap = 200, 
                 optim.method = "nlminb")

modindices(fitIGD, sort. = T)

fitIgdResults <- list(
  summary = summary(fitIGD, standardized = TRUE),
  fitMeasures = fitmeasures(fitIGD, fitMeasures),
  factorScores = parameterestimates(fitIGD, standardized = TRUE)
)

#############
igdModel <- ' 
      # latent variable
      iIgd =~ 1*igd9sf_Latent + 1*igd9sf_Latent.w2 + 1*igd9sf_Latent.w3
      sIgd =~ 0*igd9sf_Latent + 1*igd9sf_Latent.w2 + 2*igd9sf_Latent.w3
      
      iIgd ~ adhd_Latent + net_addiction_Latent + soc_media_addiction_Latent + impulsivity_Latent + stress_Latent + aggression_hostility_Latent + escape_Latent + social_support_Latent + selfassessment_Latent + anxiety_Latent + depression_Latent + social_anxiety_Latent
      sIgd ~ adhd_Latent + net_addiction_Latent + soc_media_addiction_Latent + impulsivity_Latent + stress_Latent + aggression_hostility_Latent + escape_Latent + social_support_Latent + selfassessment_Latent + anxiety_Latent + depression_Latent + social_anxiety_Latent
      
      iIgd ~~ sIgd
      
      # lag1 autocorrelations
      # igd9sf_Latent ~~ igd9sf_Latent.w2
      # igd9sf_Latent.w2 ~~ igd9sf_Latent.w3
      '

fitIGD <- growth(igdModel, data = dat, missing = "fiml", mimic = "Mplus", bootstrap = nBoot)
summary(fitIGD)

