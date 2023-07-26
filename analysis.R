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
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, fig.align = "center")
rm(list = ls())

# TODO Qs
# Zatial nevieme: Cross-Lagged Panel Analysis: How do the variables influence each other over time? For example, does impulsivity at wave 1 predict IGD symptoms at wave 2, controlling for IGD symptoms at wave 1?
# Check velke kovariancie pri b-ss networku
# lvnet package <-  tam su definovane nazvy matrixov

# Define fitmeasures
fitMea <- c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr")

# Define the number of bootstrap samples
nBoot <- 1000

# Source scripts
source("functions.R")
source("dataWrangling.R")

# Variable names
predictors <- c("adhd_Latent", "net_addiction_Latent", "soc_media_addiction_Latent", 
                "impulsivity_Latent", "stress_Latent", "aggression_hostility_Latent", 
                "escape_Latent", "social_support_Latent", "anxiety_Latent", 
                "depression_Latent", "social_anxiety_Latent")

# Create a list of symptom variable names
symptoms <- list(
  igd = paste0("igd9sf_", 1:9),
  gd = paste0("gdt_", 1:4),
  selfassessment = paste0("selfassessment_", 1:2)
)

# Define the GD latent variable names
gdVarNames <- c("igd9sf_Latent", "gdt_Latent", "selfassessment_Latent")

# Descriptives ------------------------------------------------------------
#'# Descriptives
# Calculate and prints the correlation matrix (+ Means, SDs) for each wave
waves <- list("MeanScore", "MeanScore.w2", "MeanScore.w3")
names(waves) <- paste0("wave", 1:3)

waveCorrMatrices <- lapply(waves, function(wave) {
  data <- dat %>% select(ends_with(wave))
  names(data) <- gsub(wave, "", names(data))
  getCorrelationMatrix(data)
})
names(waveCorrMatrices) <- names(waves)
waveCorrMatrices

#'## Percentages of GD cases
# Subset the data for 'igd9sf' and 'gdt' variables
igd9sfVars <- grep("igd9sf", names(dat), value = TRUE)
igd9sfVars <- igd9sfVars[!grepl("_Latent|MeanScore", igd9sfVars)]
igd9sfData <- dat[, igd9sfVars]

gdtVars <- grep("gdt", names(dat), value = TRUE)
gdtVars <- gdtVars[!grepl("_Latent|MeanScore", gdtVars)]
gdtData <- dat[, gdtVars]

selfassessment1Vars <- grep("selfassessment_1", names(dat), value = TRUE)
selfassessment1Vars <- selfassessment1Vars[!grepl("_Latent|MeanScore", selfassessment1Vars)]
selfassessment1Data <- dat[, selfassessment1Vars]

selfassessment2Vars <- grep("selfassessment_2", names(dat), value = TRUE)
selfassessment2Vars <- selfassessment2Vars[!grepl("_Latent|MeanScore", selfassessment2Vars)]
selfassessment2Data <- dat[, selfassessment2Vars]

# Function to compute the number of participants scoring >3 concurrently in all variables
countCases <- function(data) {
  wave1Vars <- grep("^[^.]*$", names(data), value = TRUE)
  wave2Vars <- grep("\\.w2$", names(data), value = TRUE)
  wave3Vars <- grep("\\.w3$", names(data), value = TRUE)
  
  wave1Scores <- rowSums(data[, wave1Vars, drop = FALSE] > 3, na.rm = TRUE) == length(wave1Vars)
  wave2Scores <- rowSums(data[, wave2Vars, drop = FALSE] > 3, na.rm = TRUE) == length(wave2Vars)
  wave3Scores <- rowSums(data[, wave3Vars, drop = FALSE] > 3, na.rm = TRUE) == length(wave3Vars)
  
  rq1 <- c(sum(wave1Scores, na.rm = TRUE)/length(wave1Scores), sum(wave2Scores, na.rm = TRUE)/length(wave2Scores), sum(wave3Scores, na.rm = TRUE) / length(wave3Scores)) * 100
  rq2 <- sum(wave1Scores | wave2Scores | wave3Scores, na.rm = TRUE) / length(wave2Scores) * 100
  rq3 <- sum(apply(cbind(wave1Scores, wave2Scores, wave3Scores), 1, combineTrueFalse), na.rm = TRUE) / length(wave3Scores) * 100
  
  list("Above cut-off per waves" = rq1, "Above cut-off in any wave" = rq2, "Above cut-off in all waves" = rq3)
}

# Apply the function to the igd9sf and gdt datasets
(gdPercentage <- list(
  igd9sf = countCases(igd9sfData),
  gdt = countCases(gdtData),
  selfassessment1 = countCases(selfassessment1Data),
  selfassessment2 = countCases(selfassessment2Data)
))

# Reliability -------------------------------------------------------------
#'# Reliability
reliabilities <- lapply(dataList, computeReliability, variableList = variableList[-c(11:14)])
reliabilities_df <- do.call(rbind, lapply(reliabilities, data.frame))

# Create a list to store the two dataframes
reliabilities_list <- list()

# Split the dataframe into two: one for alpha and one for omega
reliabilities_list$alpha <- t(reliabilities_df[, grep("cronbachAlpha", names(reliabilities_df))])
reliabilities_list$omega <- t(reliabilities_df[, grep("omega", names(reliabilities_df))])
reliabilities_list <- lapply(reliabilities_list, function(df) {
  rownames(df) <- str_replace_all(rownames(df), ".cronbachAlpha", "")
  return(df)
})

(reliabilityEstimates <- round_df(as.data.frame(reliabilities_list), 3))

# CFA results -------------------------------------------------------------
#'# Confirmatory factor analyses
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
round_df(summaryTableWide, 3) %>% select(variable, wave, chisq.scaled, df.scaled, pvalue.scaled, cfi.scaled, tli.scaled, rmsea.scaled, rmsea.ci.lower.scaled, rmsea.ci.upper.scaled, srmr)

# Longitudinal Changes ----------------------------------------------------
#'# Longitudinal changes
#'
#'## Longitudinal interplay between gaming disorder and various predictors
#'
#' Ignore the models throwing out an error.
# Function for fitting the alternative model if the primary model returns an error or warning

alternative_model <- function(gdVar, var) {
  tryCatch({
    # alternative model
    igdModel <- paste(' 
      # latent variable
      iIgd =~ 1*',gdVar,' + 1*',gdVar,'.w2 + 1*',gdVar,'.w3
      sIgd =~ 0*',gdVar,' + 1*',gdVar,'.w2 + 2*',gdVar,'.w3
      
      iPredictor =~ 1*',var,' + 1*',var,'.w2 + 1*',var,'.w3
      sPredictor =~ 0*',var,' + 1*',var,'.w2 + 2*',var,'.w3
      
      iIgd ~ iPredictor
      sIgd ~ iPredictor + sPredictor
      
      iIgd ~~ sIgd
      iPredictor ~~ sPredictor
    ', collapse = "")
    
    fitIGD <- growth(igdModel, data = dat, missing = "fiml", mimic = "Mplus", estimator = "MLR", se = "robust",
                     bootstrap = nBoot, optim.method = "nlminb")
    
    list(
      summary = summary(fitIGD, standardized = TRUE),
      fitMeasures = fitmeasures(fitIGD, fitMea),
      parameterEst = parameterestimates(fitIGD, standardized = TRUE)
    )
  }, error = function(e) {
    message(paste("Error in alternative model with gdVar: ", gdVar, " and predictor: ", var))
    message("Error message: ", e$message)
    return(NULL) # return NULL for this model if there's an error
  })
}

# Create an empty list to store results
fitResults <- list()

# Iterate over the gdVarNames
for (gdVar in gdVarNames) {
  
  tempResults <- lapply(predictors, function(var) {
    withCallingHandlers({
      igdModel <- paste(' 
        # latent variable
        iIgd =~ 1*',gdVar,' + 1*',gdVar,'.w2 + 1*',gdVar,'.w3
        sIgd =~ 0*',gdVar,' + 1*',gdVar,'.w2 + 2*',gdVar,'.w3
        
        iPredictor =~ 1*',var,' + 1*',var,'.w2 + 1*',var,'.w3
        sPredictor =~ 0*',var,' + 1*',var,'.w2 + 2*',var,'.w3
        
        iIgd ~ iPredictor
        sIgd ~ iPredictor + sPredictor
        
        iIgd ~~ sIgd
        iPredictor ~~ sPredictor
      ', collapse = "")
      
      # Fit the primary model
      fitIGD <- growth(igdModel, data = dat, missing = "fiml", mimic = "Mplus", estimator = "MLR", se = "robust",
                       bootstrap = nBoot, optim.method = "BFGS")
      
      # Return the results of the primary model
      list(
        summary = summary(fitIGD, standardized = TRUE),
        fitMeasures = fitmeasures(fitIGD, fitMea),
        parameterEst = parameterestimates(fitIGD, standardized = TRUE)
      )
    },
    # If there is an error, switch to the alternative model
    error = function(e) {
      message(paste("Error in model with gdVar: ", gdVar, " and predictor: ", var, ". Trying alternative model..."))
      alternative_model(gdVar, var)
    },
    # If there is a warning, switch to the alternative model unless the warning message is about empty cases
    warning = function(w) {
      if(!grepl("lavaan WARNING: some cases are empty and will be ignored", w$message)) {
        message(paste("Warning in model with gdVar: ", gdVar, " and predictor: ", var, ". Trying alternative model..."))
        invokeRestart("muffleWarning")
        alternative_model(gdVar, var)
      }
    })
  })
  
  # Add the results to the main results list
  names(tempResults) <- predictors
  fitResults[[gdVar]] <- tempResults
}

summaryTables <- lapply(fitResults, function(fitResult) {
  map_dfr(fitResult, function(result) {
    if (is.null(result) || !is.list(result)) {
      return(setNames(data.frame(rep(NA, 14)), c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr", "iIgd_iPred", "sIgd_iPred", "sIgd_sPred", "iIgd_iPredP", "sIgd_iPredP", "sIgd_sPredP", "iIgd_iPredP_adj")))
    }
    
    fitIndices <- result$fitMeasures
    parameters <- result$parameterEst
    
    iIgd_iPredP <- parameters$pvalue[parameters$lhs == "iIgd" & parameters$rhs == "iPredictor"]
    sIgd_iPredP <- parameters$pvalue[parameters$lhs == "sIgd" & parameters$rhs == "iPredictor"]
    sIgd_sPredP <- ifelse(length(parameters$pvalue[parameters$lhs == "sIgd" & parameters$rhs == "sPredictor"]) > 0, parameters$pvalue[parameters$lhs == "sIgd" & parameters$rhs == "sPredictor"], NA)
    
    data.frame(
      chiSq = fitIndices["chisq.scaled"],
      df = fitIndices["df.scaled"],
      pvalue = fitIndices["pvalue.scaled"],
      CFI = fitIndices["cfi.scaled"],
      TLI = fitIndices["tli.scaled"],
      RMSEA = fitIndices["rmsea.scaled"],
      SRMR = fitIndices["srmr"],
      iIgd_iPred = parameters$std.all[parameters$lhs == "iIgd" & parameters$rhs == "iPredictor"],
      sIgd_iPred = parameters$std.all[parameters$lhs == "sIgd" & parameters$rhs == "iPredictor"],
      sIgd_sPred = ifelse(length(parameters$std.all[parameters$lhs == "sIgd" & parameters$rhs == "sPredictor"]) > 0, parameters$std.all[parameters$lhs == "sIgd" & parameters$rhs == "sPredictor"], NA),
      iIgd_iPredP = iIgd_iPredP,
      sIgd_iPredP = sIgd_iPredP,
      sIgd_sPredP = sIgd_sPredP,
      iIgd_iPredP_adj = p.adjust(iIgd_iPredP, method = "holm", n = length(predictors)),
      sIgd_iPredP_adj = p.adjust(sIgd_iPredP, method = "holm", n = length(predictors)),
      sIgd_sPredP_adj = p.adjust(sIgd_sPredP, method = "holm", n = length(predictors))
    )
  }, .id = "Variable") %>% remove_rownames()
})

# Add the variable names from the original list to the summary tables
names(summaryTables) <- names(fitResults)

#' **Brief explanation:**  
#' **i prefix** represents the intercept (sometimes also referred to as the initial status or level) of the growth factor. It captures the starting point or initial level of the variable.  
#' **s prefix** represents the slope (or rate of change) of the growth factor. It captures the rate of change over time for the variable.  
#' **iIgd_iPred:** regression of the intercept of IGD (the initial level of IGD) on the intercept of the predictor (iPredictor). Essentially, this path is asking whether the initial level of the predictor is associated with the initial level of IGD.  
#' **sIgd_iPred:** regression of the slope of IGD (the rate of change in IGD) on the intercept of the predictor (iPredictor). This path is asking whether the initial level of the predictoris associated with the rate of change in IGD over time.  
#' **sIgd_sPred:** regression of the slope of IGD on the slope of the predictor (sPredictor). It is asking whether the rate of change in the predictor is associated with the rate of change in Igd over time.  
#' All parameters accompanied by their associated p-values (raw and adjusted by Holm's method).

# Results tables rounded
lapply(summaryTables, round_df, 3)

#############
#'## Gaming disorder predicted by the predictors at baseline
# Initialize the results list
t1PredResults <- list()

# Outer loop over the main variables
for (gdVar in gdVarNames) {
  
  # Inner loop over the predictors
  tempResults <- lapply(predictors, function(var) {
    igdModel <- paste(' 
        # latent variable
        iIgd =~ 1*',gdVar,' + 1*',gdVar,'.w2 + 1*',gdVar,'.w3
        sIgd =~ 0*',gdVar,' + 1*',gdVar,'.w2 + 2*',gdVar,'.w3
        
        iIgd ~ ',var,'
        sIgd ~ ',var,'
        
        iIgd ~~ sIgd
        
        # lag1 autocorrelations
        # ',gdVar,' ~~ ',gdVar,'.w2
        # ',gdVar,'.w2 ~~ ',gdVar,'.w3
      ', collapse = "")
    
    fitIGD <- growth(igdModel, data = dat, missing = "fiml", mimic = "Mplus", estimator = "MLR", se = "robust",
                     bootstrap = nBoot, optim.method = "BFGS")
    
    list(
      summary = summary(fitIGD, standardized = TRUE),
      fitMeasures = fitmeasures(fitIGD, fitMea),
      parameterEst = parameterestimates(fitIGD, standardized = TRUE)
    )
  })
  
  # Assign names to the second level of list
  names(tempResults) <- predictors
  
  # Add the results to the main results list
  t1PredResults[[gdVar]] <- tempResults
}

# Assign names to the first level of list
names(t1PredResults) <- c("igd9sf", "gdt", "selfassessment")

t1PredSummaryTableList <- lapply(t1PredResults, function(mainVarResults) {
  # Use map_dfr to create a data frame for each main variable
  mainVarTable <- purrr::map_dfr(mainVarResults, function(result) {
    if (is.null(result) || !is.list(result)) {
      return(setNames(data.frame(rep(NA, 14)), c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr", "iIgd_pred", "sIgd_pred", "iIgd_predP", "sIgd_predP", "iIgd_predP_adj", "sIgd_predP_adj")))
    }
    
    fitIndices <- result$fitMeasures
    parameters <- result$parameterEst
    
    iIgd_predP <- parameters$pvalue[parameters$lhs == "iIgd" & parameters$op == "~"]
    sIgd_predP <- parameters$pvalue[parameters$lhs == "sIgd" & parameters$op == "~"]
    
    data.frame(
      chiSq = fitIndices["chisq.scaled"],
      df = fitIndices["df.scaled"],
      pvalue = fitIndices["pvalue.scaled"],
      CFI = fitIndices["cfi.scaled"],
      TLI = fitIndices["tli.scaled"],
      RMSEA = fitIndices["rmsea.scaled"],
      SRMR = fitIndices["srmr"],
      iIgd_pred = parameters$std.all[parameters$lhs == "iIgd" & parameters$op == "~"],
      sIgd_pred = parameters$std.all[parameters$lhs == "sIgd" & parameters$op == "~"],
      iIgd_predP = iIgd_predP,
      sIgd_predP = sIgd_predP,
      iIgd_predP_adj = p.adjust(iIgd_predP, method = "holm", n = length(predictors)),   
      sIgd_predP_adj = p.adjust(sIgd_predP, method = "holm", n = length(predictors))
    )
  }, .id = "Variable") %>% remove_rownames()
  
  # Add the variable names from the original list to the summary table
  mainVarTable$Variable <- names(mainVarResults)
  
  # Return the summary table for the main variable
  mainVarTable
})

# Assign names to the list of summary tables
names(t1PredSummaryTableList) <- names(t1PredResults)

# Results tables
lapply(t1PredSummaryTableList, round_df, 3)

# Growth mixture model ----------------------------------------------------
#'# Growth mixture model
# Function to apply for each latent variable
#+ include = FALSE
gmmResults <- lapply(gdVarNames[-3], function(latentVar) {
  # Perform hlme modeling
  gmm1 <- hlme(as.formula(paste(latentVar, "~ wave")), subject = "ID", random =~ 1 + wave, ng = 1, verbose = FALSE, data = datLong %>% mutate(ID = as.numeric(as.factor(ID)))) 
  gmm2 <- gridsearch(rep = 100, maxiter = 10, minit = gmm1, hlme(as.formula(paste(latentVar, "~ wave")), subject = "ID", random =~ 1 + wave, ng = 2, verbose = FALSE, data = datLong %>% mutate(ID = as.numeric(as.factor(ID))), mixture = ~ wave, nwg = T)) 
  gmm3 <- gridsearch(rep = 100, maxiter = 10, minit = gmm1, hlme(as.formula(paste(latentVar, "~ wave")), subject = "ID", random =~ 1 + wave, ng = 3, verbose = FALSE, data = datLong %>% mutate(ID = as.numeric(as.factor(ID))), mixture = ~ wave, nwg = T)) 
  
  # Extract class memberships (posterior probabilities)
  post_probs_gmm3 <- gmm3$pprob
  
  # Get the most probable class for each subject
  classMembershipGmm3 <- post_probs_gmm3 %>% select(ID, max.col(.)) %>% setNames(c("ID", paste0(latentVar, "Class")))
  
  list(
    gmm1 = gmm1,
    gmm2 = gmm2,
    gmm3 = gmm3,
    summaryTable = summarytable(gmm1, gmm2, gmm3),
    postProb = postprob(gmm3, threshold = .5),
    estimates = estimates(gmm3),
    classMembership = classMembershipGmm3 # add class membership to the list returned
  )
})

# Assign names to the list of results
names(gmmResults) <- c("igd9sf", "gdt")

# Append the class membership variable to the data by ID
for (i in 1:length(gmmResults)) {
  dat <- left_join(dat, gmmResults[[i]]$classMembership, by = "ID")
  datLong <- left_join(datLong, gmmResults[[i]]$classMembership, by = "ID")
}

#+ include = TRUE
#'## Results of GMMs
# Extract the 'summaryTable' from each element in the list
(summaryTableList <- lapply(gmmResults, function(x) x$summaryTable))

# Extract the first element from the 'postProb' list from each element in the list
(postProbList <- lapply(gmmResults, function(x) x$postProb[[1]]))

# Extract the 'estimates' from each element in the list
(estimatesList <- lapply(gmmResults, function(x) x$estimates))

# Names of characteristics to compare across groups and names of classes
characteristics <- c("gender", "age", "gaming_time", "multiplayer")
classNames <- c("igd9sf_LatentClass", "gdt_LatentClass")

# Function to do the entire procedure for a given LatentClass variable
processLatentClass <- function(latentClassVariable) {
  # Applying the tests
  testResults <- map_df(characteristics, ~diffTest(dat, .x, latentClassVariable))
  
  # Redefine gender as numeric
  dataFiltered <- dat %>% mutate(gender = as.numeric(gender))
  
  # Split the data by the grouping variable
  dataSplit <- split(dataFiltered, dataFiltered[[latentClassVariable]])
  
  # Apply describe function to each subset
  describeList <- lapply(dataSplit, function(df) {
    sapply(characteristics, function(char) {
      descriptives <- describe(df[[char]])
      return(c(mean = descriptives$mean, sd = descriptives$sd))
    })
  })
  
  # Convert list to data frame
  describeDf <- cbind(do.call(rbind, describeList), class = rep(1:3, each = 2))
  
  # Return the test results and the data frame of descriptive statistics
  list(testResults = testResults, describeDf = describeDf)
}

# Apply the function to both "igd9sf_LatentClass" and "gdt_LatentClass"
resultsGMM <- lapply(classNames, processLatentClass)
names(resultsGMM) <- classNames
resultsGMM

#'## Plot the trendlines for the distinct subgroups differing in their GD trajectory
# Find the variable names in the dataset that contain the string "_LatentClass"
latentClassVarNames <- names(datLong)[grepl("_LatentClass", names(datLong))]

# Function to create a plot for a given GD variable and a given LatentClass variable
plotTrends <- function(gdVariable, latentClassVariable) {
  data <- datLong %>% filter(!is.na(datLong[[latentClassVariable]]))
  varName <- toupper(gsub("_Latent", "", gdVariable))
  latentClassName <- toupper(gsub("_LatentClass", "", latentClassVariable))
  
  # Get the percentages for the current LatentClass variable from postProbList
  percentages <- round(postProbList[[gsub("_LatentClass", "", latentClassVariable)]][2, ], 0)
  # Create a named vector of labels for the factor levels
  factorLabels <- paste0(1:3, " (", percentages, "%)")
  names(factorLabels) <- 1:3
  
  ggplot(data, aes_string(x = "wave", y = gdVariable, group = latentClassVariable)) +
    geom_smooth(aes(color = as.factor(get(latentClassVariable))), method = "loess", span = .4, se = FALSE, size = 1.5) +
    labs(title = paste("Classification based on", latentClassName),
         x = "Wave",
         y = varName,
         color = "Latent class") +
    scale_x_continuous(breaks = 1:3) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    scale_color_discrete(labels = factorLabels) +  # Set the labels in the legend
    theme_minimal()
}

# Create the plots
gmmPlots <- list()
for (gdVar in gdVarNames[-3]) {
  for (lcVar in latentClassVarNames) {
    gmmPlots[[paste(gdVar, lcVar, sep = "_")]] <- plotTrends(gdVar, lcVar)
  }
}

# Arrange the plots in a grid
grid.arrange(gmmPlots[[1]], gmmPlots[[2]], gmmPlots[[3]], gmmPlots[[4]], ncol = 2)

# Development of GD symptoms over time ---------------------------------------
#'# Development of GD symptoms over time
#'
#+ include = FALSE
lmerModels <- function(symptoms, datLong) {
  resultsSymptomsList <- list()
  
  for (i in seq_along(symptoms)) {
    
    resultsSymptoms <- list()
    
    for (symptom in symptoms[[i]]) {
      null_formula <- as.formula(paste0(symptom, " ~ 1 + (1 + wave|ID)"))
      full_formula <- as.formula(paste0(symptom, " ~ wave + (1 + wave|ID)"))
      
      null_model <- lmer(null_formula, data = datLong, REML = T)
      full_model <- lmer(full_formula, data = datLong, REML = T)
      
      lr_test <- anova(null_model, full_model)
      
      wave_coef <- fixef(full_model)["wave"]
      wave_se <- sqrt(diag(vcov(full_model))["wave"])
      random_effects <- as.data.frame(VarCorr(full_model)$ID)
      
      intercept_sd <- sqrt(random_effects[1, 1])
      wave_sd <- sqrt(random_effects[2, 2])
      residual_sd <- attr(VarCorr(full_model), "sc")
      intercept_wave_corr <- random_effects[2, 1] / (intercept_sd * wave_sd)
      
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
      
      symptom_list <- list(
        null_model = null_model,
        full_model = full_model,
        null_summary = summary(null_model),
        full_summary = summary(full_model),
        symptom_results = symptom_results
      )
      
      resultsSymptoms[[symptom]] <- symptom_list
    }
    
    summaryTableSymptoms <- do.call(rbind, lapply(resultsSymptoms, function(x) x$symptom_results)) %>% remove_rownames()
    names(summaryTableSymptoms) <- c("Symptom", "Wave Coefficient", "SE", "Intercept SD", "Wave SD", "Residual SD", "Correlation", "p-value")
    summaryTableSymptoms$Symptom <- gsub("igd9sf_", "IGDS9-SF ", summaryTableSymptoms$Symptom)
    
    resultsSymptomsList[[i]] <- list(
      results = resultsSymptoms,
      summaryTable = summaryTableSymptoms
    )
  }
  
  names(resultsSymptomsList) <- c("igd9sf", "gdt", "selfassessment")
  
  return(resultsSymptomsList)
}

resultsSymptomsList <- lmerModels(symptoms, datLong)

#+ include = TRUE
# Extract the 'summaryTable' from each element in the list
#'## Summary tables for the lmer models
#'
#' **Wave Coefficient:** the estimated fixed effect of wave on the symptom, indicating the average change in the symptom per unit increase in wave.  
#' **Standard Error:** the standard error of the wave coefficient, providing a measure of its uncertainty.  
#' **Intercept SD and Wave SD:** the standard deviations of the random intercepts and slopes, respectively, indicating the variability in these effects across participants.  
#' **Intercept-Wave Correlation:** the correlation between the random intercepts and slopes, indicating the extent to which participants with higher baseline symptom levels also tend to show larger wave effects.  
#' **LR Test p-value:** the p-value from the likelihood ratio test comparing the full model (including wave as a fixed effect and random slope) to the null model (excluding wave), providing evidence for whether wave has a significant effect on the symptom.  

(symptomsSummaryTableList <- lapply(resultsSymptomsList, function(x) round_df(x$summaryTable, 3)))
#' The negative correlations suggest that individuals who start with higher symptom levels tend to experience a smaller change or improvement in their symptoms over time.
#' 
#'### Mean wave coefficient for igd9sf
median(symptomsSummaryTableList$igd9sf$`Wave Coefficient`)
#'### Mean ratio of intercept to wave SDs for igd9sf
median(symptomsSummaryTableList$igd9sf$`Intercept SD`/symptomsSummaryTableList$igd9sf$`Wave SD`)

#'#### Plot
# Plot the fitted values
# Extract the fitted values from each model
# Apply a function to each element of the 'symptoms' list
gdPlots <- lapply(seq_along(symptoms), function(i) {
  
  fitted_values <- data.frame()
  
  for (symptom in symptoms[[i]]) {
    # Get the data used in the model
    model_data <- getME(resultsSymptomsList[[i]]$results[[symptom]][["full_model"]], "X")
    
    # Extract the fitted values
    symptom_fitted_values <- data.frame(
      wave = as.data.frame(model_data)$wave,
      symptom = rep(gsub("igd9sf_", "IGDS9-SF ", gsub("gdt_", "GDT ", symptom)), nrow(model_data)),
      fitted_value = fitted(resultsSymptomsList[[i]]$results[[symptom]][["full_model"]])
    )
    
    # Add the fitted values to the data frame
    fitted_values <- rbind(fitted_values, symptom_fitted_values)
  }
  
  # Calculate the mean fitted value for each wave and symptom
  mean_fitted_values <- aggregate(fitted_value ~ wave + symptom, data = fitted_values, FUN = mean)
  
  # Create the plot
  plot <- ggplot(mean_fitted_values, aes(x = wave, y = fitted_value, color = symptom)) +
    geom_smooth(method = "loess", se = FALSE, span = .7, size = 1.5) +
    scale_x_continuous(breaks = c(1, 2, 3)) +
    coord_cartesian(ylim = c(1, 3.2)) + 
    scale_color_brewer(palette = "Set1") +
    labs(x = "Wave", y = "Mean fitted value", color = "Symptom") +
    theme_minimal()
  
  return(plot)
})

# Arrange the plots in a grid
grid.arrange(gdPlots[[1]], gdPlots[[2]], ncol = 2)

# Development of risk factors over time -----------------------------------
#'# Development of risk factors over time
# Initialize an empty list to store the results

predictorsMeans <- c("adhdMeanScore", "netAddictionMeanScore", "socMediaAddictionMeanScore", 
                "impulsivityMeanScore", "stressMeanScore", "aggressionHostilityMeanScore", 
                "escapeMeanScore", "socialSupportMeanScore", "anxietyMeanScore", 
                "depressionMeanScore", "socialAnxietyMeanScore")

resultsPredictors <- list()

# Fit separate mixed-effects models for each risk factor
for (predictor in predictorsMeans) {
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

# Print the table 
#'## Summary tables for the lmer models
round_df(summaryTablePredictors, 3)

#'#### Plot
# Plot the fitted values
# Extract the fitted values from each model
predFitted <- data.frame()
for (predictor in predictorsMeans) {
  # Get the data used in the model
  model_data <- getME(resultsPredictors[[predictor]][["full_model"]], "X")

  # Extract the fitted values
  predictorFittedValues <- data.frame(
    wave = as.data.frame(model_data)$wave,
    predictor = rep(predictor, nrow(model_data)),
    fitted_value = fitted(resultsPredictors[[predictor]][["full_model"]])
  )

  # Add the fitted values to the data frame
  predFitted <- rbind(predFitted, predictorFittedValues)
}

# Calculate the mean fitted value for each wave and symptom
predMeanFitted <- aggregate(fitted_value ~ predictor + wave, data = predFitted, FUN = mean)

# Create a named vector to map predictors to their names
predictorsMeansNames <- c("ADHD", "Internet addiction", "Social media addiction", "Impulsivity",
                          "Stress", "Aggression/hostility", "Escape", "Social support", "Anxiety",
                          "Depression", "Social anxiety")
nameMap <- setNames(predictorsMeansNames, predictorsMeans)

# Replace predictor names with their corresponding names in the data
predMeanFitted$predictor <- nameMap[predMeanFitted$predictor]

# Create the plot
ggplot(predMeanFitted, aes(x = wave, y = fitted_value, color = predictor)) +
  geom_smooth(method = "loess", se = FALSE, span = .8, size = 1) +
  scale_x_discrete(breaks = c("1", "2", "3"), expand = c(0,0)) +
  labs(x = "Wave", y = "Mean fitted value", color = "Risk factor") +
  theme_minimal()

# Panel VAR network model -------------------------------------------------
#'# Panel VAR network model
#'
#' **For IGD9-SF and GDT, respectively.**
alpha <- 0.01
adjust <- "none"
searchstrategy <- "stepup"

# List of symptom prefixes
symptom_prefixes <- c("igd9sf", "gdt")
psychonetricsResults <- list()

# Function to run the model for a certain symptom prefix
psychonetricsEstimation <- function(symptom_prefix) {
  # Select data columns that start with the symptom prefix
  datVAR <- dat %>% select(starts_with(symptom_prefix) & !contains("Latent") & !contains("MeanScore"))
  nVAR <- nrow(datVAR)
  covVAR <- ((nVAR - 1)/nVAR) * cov(datVAR, use = "pairwise.complete.obs")
  meansVAR <- colMeans(datVAR, na.rm = TRUE)
  
  # Define design
  designVAR <- matrix(dat %>% select(starts_with(symptom_prefix) & !contains("Latent")  & !contains("MeanScore")) %>% colnames(),
                      nrow = ncol(dat %>% select(starts_with(symptom_prefix) & !contains("Latent")  & !contains("MeanScore")))/3,
                      ncol = 3)
  
  # Run the panelgvar model
  model <- panelgvar(covs = covVAR, nobs = nVAR, vars = designVAR, estimator = "ML") %>%
    runmodel(bounded = FALSE)
  
  modelPruned <- model %>%
    runmodel %>%
    prune(alpha = alpha, adjust = adjust, recursive = TRUE)
  
  # Search strategy
  if (searchstrategy == "stepup"){
    modelPruned <- modelPruned %>%  stepup(alpha = alpha, criterion = "bic")
  } else if (searchstrategy == "modelsearch"){
    modelPruned <- modelPruned %>%  modelsearch(prunealpha = alpha, addalpha = alpha)
  }
  
  model %>% fit
  modelPruned %>% fit
  
  # Comparison
  (comp <- compare(original = model, pruned = modelPruned))
  
  # Differences:
  comp$AIC[1] - comp$AIC[2]
  comp$BIC[1] - comp$BIC[2]
  
  # Model fit
  fitModel <- model %>% fit
  fitModelPruned <- modelPruned %>% fit
  model_params <- model %>% parameters
  modelPruned_params <- modelPruned %>% parameters
  
  # Labels
  labels <- paste0(symptom_prefix, seq_along(dat %>% select(starts_with(symptom_prefix) & !contains("Latent") & !contains(".w") & !contains("MeanScore")) %>% colnames()))
  
  # Get matrices
  (temporal <- temporalCov <- getmatrix(modelPruned, "PDC"))
  (contemporaneous <- contemporaneousCov <- getmatrix(modelPruned, "omega_zeta_within"))
  (between <- betweenCov <- getmatrix(modelPruned, "omega_zeta_between"))

  # (between <- betweenCov <- getmatrix(modelPruned, "lowertri_epsilon_between"))
  
  # Plot networks
  layout(t(1:3))
  max <- max(c(abs(temporal), abs(contemporaneous), abs(between)))
  qgraph(temporal, layout = "circle", labels = toupper(labels), directed = TRUE, minimum = .1, maximum = max, mar = rep(5,4),
         title = "Temporal network", shape = "rectangle", theme = "colorblind", label.scale.equal = TRUE,
         asize = 8, vsize = 13, label.cex = 1.5)
  box("figure")
  
  qgraph(contemporaneous, layout = "circle", labels = toupper(labels), minimum = .1, maximum = max, mar = rep(5,4),
         title = "Contemporaneous network", shape = "rectangle", theme = "colorblind", label.scale.equal = TRUE,
         asize = 8, vsize = 13, label.cex = 1.5)
  box("figure")
  
  qgraph(between, layout = "circle", labels = toupper(labels), minimum = .1, maximum = max, mar = rep(5,4),
         title = "Between-subjects network", shape = "rectangle", theme = "colorblind", label.scale.equal = TRUE,
         asize = 8, vsize = 13, label.cex = 1.5)
  box("figure")
  
  # Parameter estimates
  rownames(temporalCov) <- colnames(temporalCov) <- rownames(contemporaneousCov) <- colnames(contemporaneousCov) <- rownames(betweenCov) <- colnames(betweenCov) <- labels
  # Temporal
  temporalCor <- cov2cor(getmatrix(modelPruned, "PDC"))
  temporalCov[lower.tri(temporalCov)] <- temporalCor[lower.tri(temporalCor)]
  temporalMat <- round_df(as.data.frame(temporalCov), 3)
  # Contemporaneous
  contemporaneousCor <- cov2cor(getmatrix(modelPruned, "omega_zeta_within"))
  contemporaneousCov[lower.tri(contemporaneousCov)] <- contemporaneousCor[lower.tri(contemporaneousCor)]
  contemporaneousMat <- round_df(as.data.frame(contemporaneousCov), 3)
  # Between
  betweenCor <- cov2cor(getmatrix(modelPruned, "omega_zeta_between"))
  betweenCov[lower.tri(betweenCov)] <- betweenCor[lower.tri(betweenCor)]
  betweenMat <- round_df(as.data.frame(betweenCov), 3)

  # Save the results
  psychonetricsResults[[symptom_prefix]] <<- list(
    model = model,
    modelPruned = modelPruned,
    comp = comp,
    fitModel = fitModel,
    fitModelPruned = fitModelPruned,
    modelParams = model_params,
    modelPrunedParams = modelPruned_params,
    temporal = temporal,
    contemporaneous = contemporaneous,
    between = between
  )
}
# Run the function for each symptom prefix
psychonetricsModels <- lapply(symptom_prefixes, psychonetricsEstimation)
