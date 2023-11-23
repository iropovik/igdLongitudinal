# rm(list = ls())

list.of.packages <- c("readxl", "lavaan", "vars", "psychonetrics", "qgraph", "lme4", "lcmm", "tidyverse", "magrittr", "psych", "plyr", "semTools", "GPArotation", "careless", "knitr", "kableExtra", "gplots", "gridExtra", "naniar", "broom", "ggplot2", "gridExtra", "parallel")
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
if (!exists("fitMea")) {
  fitMea <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
}

# Load data
dataList <- list(
  data1 = read_excel("data_1.xlsx"),
  data2 = read_excel("data_2.xlsx"),
  data3 = read_excel("data_3.xlsx")
)

# Remove completely empty rows from each dataset in dataList
# dataList <- lapply(dataList, function(df) df[rowSums(is.na(df)) != ncol(df), ])

# Create lookup_table and add numeric_ID to each dataset in dataList
lookup_table <- unique(unlist(sapply(dataList, `[[`, "ID")))
dataList <- lapply(dataList, function(df) transform(df, ID = match(ID, lookup_table)))

# Reverse-scaling of items
reverseItems <- c("stress_4", "stress_5", "stress_7", "stress_8")
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
  
  # Remove participants who missed both attention checks or 
  # missed one of the attention checks and was a significant multivariate outlier
  data <- data %>% filter((control_item_1 == 5 & control_item_2 == 0) | 
                            (control_item_1 != 5 & mahalanobis > cutoff) |
                            (control_item_2 != 0 & mahalanobis > cutoff))
})

# Compute mean scores
dataList <- lapply(dataList, function(data) {
  data <- data %>% mutate(
    igd9sfMeanScore = rowMeans(select(data, igd9sf_1:igd9sf_9), na.rm = TRUE),
    gdtMeanScore = rowMeans(select(data, gdt_1:gdt_4), na.rm = TRUE),
    selfAssessmentMeanScore = rowMeans(select(data, selfassessment_1:selfassessment_2), na.rm = TRUE),
    adhdMeanScore = rowMeans(select(data, adhd_1:adhd_6), na.rm = TRUE),
    netAddictionMeanScore = rowMeans(select(data, net_addiction_1:net_addiction_6), na.rm = TRUE),
    socMediaAddictionMeanScore = rowMeans(select(data, soc_media_addiction_1:soc_media_addiction_6), na.rm = TRUE),
    impulsivityMeanScore = rowMeans(select(data, impulsivity_1:impulsivity_10), na.rm = TRUE),
    stressMeanScore = rowMeans(select(data, stress_1:stress_10), na.rm = TRUE),
    aggressionHostilityMeanScore = rowMeans(select(data, aggression_hostility_1:aggression_hostility_6), na.rm = TRUE),
    escapeMeanScore = rowMeans(select(data, escape_1:escape_4), na.rm = TRUE),
    socialSupportMeanScore = rowMeans(select(data, social_support_1:social_support_6), na.rm = TRUE),
    anxietyMeanScore = rowMeans(select(data, anxiety_1:anxiety_2), na.rm = TRUE),
    depressionMeanScore = rowMeans(select(data, depression_1:depression_2), na.rm = TRUE),
    socialAnxietyMeanScore = rowMeans(select(data, social_anxiety_1:social_anxiety_3), na.rm = TRUE),
    selfesteemMeanScore = selfesteem,
    educationMeanScore = education
  )
  return(data)
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
  socialAnxiety = c("social_anxiety_1", "social_anxiety_2", "social_anxiety_3"),
  selfesteem = "selfesteem",
  education = "education"
)

# For scales within each dataset in the list, compute factor or principal component score
results <- lapply(dataList, cfaFunction, variableList = variableList, boot = nBoot, fitMeasures = fitMea)

# Extract the data and results lists
dataList <- lapply(results, `[[`, "data")
resultsList <- lapply(results, `[[`, "modelResults")

# Create data objects
# Wide dataset
dat <- full_join(dataList$data1, dataList$data2, by = "ID", suffix = c("", ".w2")) %>%
  full_join(dataList$data3, by = "ID", suffix = c("", ".w3"))

# Long data
datLong <- bind_rows(dataList) %>% arrange(ID, wave)

# Define gender as factor and define the categories
dat <- dat %>% 
  mutate(gender = as.factor(gender)) %>% 
  filter(gender %in% c(1, 2))

# Delete further unused variables
delvars <- names(dat) %in% c(
  "wave", "wave.w2", "wave.w3"
)
dat <- dat[!delvars]
rm(delvars)
