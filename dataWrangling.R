# rm(list = ls())
# to do gpt research questions

list.of.packages <- c("readxl", "lavaan", "vars", "psychonetrics", "qgraph", "lme4", "semPlot", "tidyverse", "magrittr", "psych", "plyr", "semTools", "careless", "knitr", "kableExtra", "gplots", "naniar", "ggplot2", "gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load required libraries
#+ include = FALSE
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)

# Define the number of bootstrap samples
if (!exists("nBoot")) {
  nBoot <- 200
}

# Define fitmeasures
if (!exists("fitMeasures")) {
  fitMeasures <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
}

# Load data
dataList <- list(
  data1 = read_excel("data_1.xlsx"),
  data2 = read_excel("data_2.xlsx"),
  data3 = read_excel("data_3.xlsx")
)

# Reverse-scaling of items
reverseItems <- c("stress_4", "stress_5", "stress_7", "stress_8")
reverseScale <- function(dataset) {
  for (item in reverseItems) {
    if (item %in% names(dataset)) {
      dataset[[item]] <- 4 - dataset[[item]]
    }
  }
  return(dataset)
}

# Apply the function to each dataset in the list
dataList <- lapply(dataList, reverseScale)

# Remove participants with gaming time == 0
dataList <- lapply(dataList, function(data) {
  data <- data %>% subset(gaming_time != 0 | is.na(gaming_time))
  return(data)
})

# Check missing data patterns
(naProportions <- lapply(dataList, function(data) {
  selectedData <- select(data, gender:social_support_6)
  naCount <- sum(is.na(selectedData))
  proportionNA <- naCount / (nrow(selectedData) * ncol(selectedData))
}))

(missingVarSummary <- lapply(dataList, miss_var_summary))
(missingCaseSummary <- lapply(dataList, miss_case_summary))

# Number of participants missing both attention checks
lapply(dataList, function(data) {
  nrow(data) - data %>% filter((control_item_1 == 5 | control_item_2 == 0)) %>% nrow()
})

# Mahalanobis distance
dataList <- lapply(dataList, function(data) {
  # The alpha quantile of the chi-square distribution was used (Rousseeuw & Van Zomeren, 1990)
  # with the alpha threshold set to 0.025 (Cabana, 2019).
  outliers <- data %>% select(age:social_support_6)
  cutoff <- qchisq(p = 1-0.025 , df = ncol(outliers))
  data$mahalanobis <- mahad(select(data, age:social_support_6),
                            confidence = .99, plot = FALSE)
  
  # Remove participants who missed both attention checks
  data <- data %>% filter((control_item_1 == 5 | control_item_2 == 0))
})

# Create factor scores or mean z-scores for variables that are indicated by less than 4 items
# Define the list of variables for each factor
variableList <- list(
  gdt = c("gdt_1", "gdt_2", "gdt_3", "gdt_4"),
  igd9sf = c("igd9sf_1", "igd9sf_2", "igd9sf_3", "igd9sf_4", "igd9sf_5", "igd9sf_6", "igd9sf_7", "igd9sf_8", "igd9sf_9"),
  adhd = c("adhd_1", "adhd_2", "adhd_3", "adhd_4", "adhd_5", "adhd_6"),
  netAddiction = c("net_addiction_1", "net_addiction_2", "net_addiction_3", "net_addiction_4", "net_addiction_5", "net_addiction_6"),
  socMediaAddiction = c("soc_media_addiction_1", "soc_media_addiction_2", "soc_media_addiction_3", "soc_media_addiction_4", "soc_media_addiction_5", "soc_media_addiction_6"),
  impulsivity = c("impulsivity_1", "impulsivity_2", "impulsivity_3", "impulsivity_4", "impulsivity_5", "impulsivity_6", "impulsivity_7", "impulsivity_8", "impulsivity_9", "impulsivity_10"),
  stress = c("stress_1", "stress_2", "stress_3", "stress_4", "stress_5", "stress_6", "stress_7", "stress_8", "stress_9", "stress_10"),
  aggressionHostility = c("aggression_hostility_1", "aggression_hostility_2", "aggression_hostility_3", "aggression_hostility_4", "aggression_hostility_5", "aggression_hostility_6"),
  escape = c("escape_1", "escape_2", "escape_3", "escape_4"),
  socialSupport = c("social_support_1", "social_support_2", "social_support_3", "social_support_4", "social_support_5", "social_support_6"),
  selfAssessment = c("selfassessment_1", "selfassessment_2"),
  anxiety = c("anxiety_1", "anxiety_2"),
  depression = c("depression_1", "depression_2"),
  socialAnxiety = c("social_anxiety_1", "social_anxiety_2", "social_anxiety_3")
)

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
                 std.lv = TRUE, mimic = "Mplus", estimator = "MLR", missing = "fiml",
                 bootstrap = nBoot
      )
      
      factorScores <- as.numeric(lavPredict(fit))
      factorScoresData <- cbind(factorScoresData, setNames(data.frame(factorScores), paste0(variableName, "_Latent")))
      
      modelResults[[variableName]] <- list(
        summary = summary(fit),
        fitMeasures = fitmeasures(fit, fitMeasures),
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

# Apply the function to each dataset in the list
results <- lapply(dataList, cfaFunction, variableList = variableList, boot = nBoot, fitMeasures = fitMeasures)

# Extract the data and results lists
dataList <- lapply(results, `[[`, "data")
resultsList <- lapply(results, `[[`, "modelResults")

# Create data objects
# Wide dataset
dat <- full_join(dataList$data1, dataList$data2, by = "ID", suffix = c("", ".w2")) %>%
  full_join(dataList$data3, by = "ID", suffix = c("", ".w3"))

# Long data
datLong <- bind_rows(dataList) %>% arrange(ID, wave)

# Delete further unused variables
delvars <- names(dat) %in% c(
  
)
dat <- dat[!delvars]
rm(delvars)

##############################
# Compute mean scores -----------------------------------------------------
# dataList <- lapply(dataList, function(data) {
#   data <- data %>% mutate(
#     gdt = rowMeans(select(data, gdt_1:gdt_4), na.rm = TRUE),
#     igd9sf = rowMeans(select(data, igd9sf_1:igd9sf_9), na.rm = TRUE),
#     selfAssessment = rowMeans(select(data, selfassessment_1:selfassessment_2), na.rm = TRUE),
#     anxiety = rowMeans(select(data, anxiety_1:anxiety_2), na.rm = TRUE),
#     depression = rowMeans(select(data, depression_1:depression_2), na.rm = TRUE),
#     socialAnxiety = rowMeans(select(data, social_anxiety_1:social_anxiety_3), na.rm = TRUE),
#     adhd = rowMeans(select(data, adhd_1:adhd_6), na.rm = TRUE),
#     netAddiction = rowMeans(select(data, net_addiction_1:net_addiction_6), na.rm = TRUE),
#     socMediaAddiction = rowMeans(select(data, soc_media_addiction_1:soc_media_addiction_6), na.rm = TRUE),
#     impulsivity = rowMeans(select(data, impulsivity_1:impulsivity_10), na.rm = TRUE),
#     stress = rowMeans(select(data, stress_1:stress_10), na.rm = TRUE),
#     aggressionHostility = rowMeans(select(data, aggression_hostility_1:aggression_hostility_6), na.rm = TRUE),
#     escape = rowMeans(select(data, escape_1:escape_4), na.rm = TRUE),
#     socialSupport = rowMeans(select(data, social_support_1:social_support_6), na.rm = TRUE)
#   )
#   return(data)
# })

# gdtModel <- 'gdtF =~ gdt_1 + gdt_2 + gdt_3 + gdt_4'
# gdtFit <- cfa(gdtModel, data,
#               std.lv = TRUE, mimic = "Mplus", estimator = "WLSMV",
#               bootstrap = boot, missing = "pairwise",
#               ordered = c("gdt_1", "gdt_2", "gdt_3", "gdt_4"))
# fitmeasures(gdtFit, fitMeasures)
# data$gdtLatent <- lavPredict(gdtFit)