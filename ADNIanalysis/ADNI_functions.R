
calc_rate_of_change <- function(sub_data, val_col, time_col) {
  
  # Ensure the value and time columns are numeric
  sub_data[[val_col]] <- as.numeric(sub_data[[val_col]])
  sub_data[[time_col]] <- as.numeric(sub_data[[time_col]])
  
  # Remove rows with NA values in the relevant columns
  sub_data <- sub_data %>%
    filter(!is.na(!!sym(val_col)) & !is.na(!!sym(time_col)))
  
  # Create the formula dynamically using the column names
  formula <- as.formula(paste(val_col, "~", time_col))
  
  # Perform linear regression using the dynamically created formula
  if (nrow(sub_data) > 1) {  # Ensure there is enough data for regression
    model <- lm(formula, data = sub_data)
    return(coef(model)[[time_col]])  # Return the slope coefficient
  } else {
    return(NA)  # Return NA if there isn't enough data
  }
}
get_RofC_df <- function(df, id_name, value_name, time_name, flip_sign = FALSE) {
  
  # Define column symbols
  subject_col <- sym(id_name)
  val_col <- value_name
  time_col <- time_name
  
  # Apply the function to each subject and get the rate of change
  RofC <- df %>%
    group_by(!!subject_col) %>%
    group_modify(~ tibble(rate_of_change = calc_rate_of_change(.x, val_col, time_col)))
  
  if (flip_sign) {
    RofC$rate_of_change <- RofC$rate_of_change * -1
  }
  
  # Compute Z-scores and remove outliers
  RofC$Z <- (RofC$rate_of_change - mean(RofC$rate_of_change, na.rm = TRUE)) / 
    sd(RofC$rate_of_change, na.rm = TRUE)
  RofC <- data.frame(RofC)
  RofC <- RofC[RofC$Z < 5 & RofC$Z > -5 & !is.na(RofC$Z), ]
  RofC <- RofC[, !names(RofC) %in% c("Z")]
  # Fit clustering model
  
  fit <- Mclust(RofC[!is.na(RofC$rate_of_change), ]$rate_of_change, G = 2, 
                model = "V")
  
  # Assign classifications
  RofC$classification <- NA
  RofC[!is.na(RofC$rate_of_change), ]$classification <- fit$classification
  
  # --- Reclassification of the left-hand tail (as in your original code) ---
  row_index <- which.min(abs(RofC$rate_of_change - median(RofC$rate_of_change, na.rm = TRUE)))
  RofC[RofC$rate_of_change < median(RofC$rate_of_change), ]$classification <- 
    RofC[row_index, ]$classification
  
  # --- Ensure class 1 is lower distribution and class 2 is upper distribution ---
  class_means <- tapply(RofC$rate_of_change, RofC$classification, mean, na.rm = TRUE)
  
  if (length(class_means) == 2 && class_means[1] > class_means[2]) {
    # If class 1 is actually the higher group, swap labels
    RofC$classification <- ifelse(RofC$classification == 1, 2, 1)
  }
  
  RofC$classification <- factor(RofC$classification, levels = c(1, 2))
  
  return(RofC)
}