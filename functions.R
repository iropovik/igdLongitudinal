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

# Function to compute reliability estimates
computeReliability <- function(data, variableList) {
  reliabilityResults <- list()
  
  for (variable in variableList) {
    # Extract the base name of the variable (remove the "_1")
    variableName <- strsplit(variable[1], "_1")[[1]][1]
    
    # Compute Cronbach's alpha and omega
    reliability <- psych::alpha(data[variable])
    cronbachAlpha <- reliability$total$raw_alpha
    omega <- psych::omega(data[variable], plot = FALSE)$omega.tot
    
    reliabilityResults[[variableName]] <- list(
      cronbachAlpha = cronbachAlpha,
      omega = omega
    )
  }
  
  return(reliabilityResults)
}
