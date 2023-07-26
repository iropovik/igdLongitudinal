cfaFunction <- function(data, variableList, boot, fitMeasures) {
  modelResults <- list()
  factorScoresData <- data.frame(ID = data$ID)
  
  for (variable in variableList) {
    # Extract the base name of the variable (remove the "_1")
    variableName <- strsplit(variable[1], "_1")[[1]][1]
    
    # Check if the column already exists in the factorScoresData data frame
    if (paste0(variableName, "_Latent") %in% names(factorScoresData)) {
      next
    }
    
    if (length(variable) >= 4) {
      # If the variable is comprised of 4 items or more, perform CFA
      modelFormula <- paste0(variable[1], "F =~ ", paste(variable, collapse = " + "))
      
      fit <- cfa(modelFormula, data,
                 std.lv = TRUE, mimic = "Mplus", estimator = "WLSMV", test = "Satterthwaite", se = "robust",
      )
      
      # Check for missing data and add NAs for these cases when calculating the factor scores
      factorScores <- rep(NA, nrow(data))
      complete_cases <- complete.cases(data[variable])
      factorScores[complete_cases] <- as.numeric(lavPredict(fit))
      
      factorScoresData <- cbind(factorScoresData, setNames(data.frame(factorScores), paste0(variableName, "_Latent")))
      
      modelResults[[variableName]] <- list(
        summary = summary(fit, standardized = TRUE),
        fitMeasures = fitmeasures(fit, fitMea),
        factorScores = factorScores
      )
    } else {
      # If the variable is comprised of less than 4 items, compute PCA score
      pcaResult <- principal(data[variable], nfactors = 1, scores = TRUE)
      pcaScore <- as.numeric(pcaResult$scores)
      factorScoresData <- cbind(factorScoresData, setNames(data.frame(pcaScore), paste0(variableName, "_Latent")))
      
      modelResults[[variableName]] <- list(
        summary = NA,
        fitMeasures = NA,
        factorScores = pcaScore
      )
    }
  }
  
  # Join the original data with the factor scores data
  data <- left_join(data, factorScoresData, by = "ID")
  return(list(data = data, modelResults = modelResults))
}

# Round data.frame function
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}

# Function to remove NA elements recursively from a nested list
removeNAElements <- function(x) {
  if (is.list(x)) {
    # Recursively apply the removeNAElements function to each element of the list
    x <- lapply(x, removeNAElements)
    
    # Remove elements with NA at the current level
    x <- x[!sapply(x, function(y) all(is.na(y)))]
  }
  return(x)
}

# Reverse-scaling of items
reverseScale <- function(dataset) {
  for (item in reverseItems) {
    if (item %in% names(dataset)) {
      dataset[[item]] <- 4 - dataset[[item]]
    }
  }
  return(dataset)
}

# Function to calculate the correlation matrix and append mean and SD
getCorrelationMatrix <- function(data) {
  corrMatrix <- corr.test(data)$r
  corrMatrix <- cbind(corrMatrix, Mean = colMeans(data, na.rm = TRUE), SD = apply(data, 2, sd, na.rm = TRUE))
  corrMatrix <- round_df(data.frame(corrMatrix), 3)  # Round the values in the matrix to 3 decimals
  
  # Change the row names
  rownames(corrMatrix) <- paste0("[", 1:nrow(corrMatrix), "] ", rownames(corrMatrix))
  # Change the column names
  colnames(corrMatrix) <- c(paste0("[", 1:(ncol(corrMatrix) - 2), "]"), "Mean", "SD")
  # Keep only the lower triangular part of the matrix, including the diagonal
  corrMatrix[, 1:(ncol(corrMatrix) - 2)][upper.tri(corrMatrix[, 1:(ncol(corrMatrix) - 2)], diag = FALSE)] <- ""
  return(corrMatrix)
}

# Function to compute reliability estimates
computeReliability <- function(data, variableList) {
  reliabilityResults <- list()
  
  for (variable in variableList) {
    # Extract the base name of the variable (remove the "_1")
    variableName <- strsplit(variable[1], "_1")[[1]][1]
    
    # Compute Cronbach's alpha and omega
    reliability <- tryCatch({
      suppressMessages(psych::omega(data[variable], plot = FALSE, nfactors = 1, poly = TRUE))
    }, error = function(e) {
      list(alpha = NA, omega.tot = NA)
    })
    
    cronbachAlpha <- reliability$alpha
    omega <- reliability$omega.tot
    
    reliabilityResults[[variableName]] <- list(
      cronbachAlpha = cronbachAlpha,
      omega = omega
    )
  }
  
  return(reliabilityResults)
}

# Function to conduct ANOVA or chi-squared test depending on the type of the variable
diffTest <- function(data, var, class) {
  if(is.numeric(data[[var]])) {
    testResult <- aov(data[[var]] ~ data[[class]], data = data) %>% tidy() %>% 
      select(term, df, statistic, p.value) %>% mutate(variable = var)
    testResult$method <- "ANOVA"
  } else {
    testResult <- chisq.test(data[[var]], data[[class]]) %>% tidy() %>% 
      select(statistic, p.value) %>% mutate(variable = var)
    testResult$method <- "Chi-squared"
  }
  return(testResult)
}

# Define a custom function to handle NA values in logical vectors for GD percentages: RQ3
combineTrueFalse <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(all(x[!is.na(x)], na.rm = TRUE))
  }
}
