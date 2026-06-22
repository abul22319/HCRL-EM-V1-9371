# monoexp_model_HCRL.R
# NOTE: this side of the script is doing the modeling itself. 
#If any lab members want changes to the calculations, they must be done here

#Fit a mono exponential model to a time series

#' @param data      Data frame containing a "Time"
#' @param variable   Name (string) of the column in `data` to fit.
#' @param direction  1 = Rise 2 = Decay
#' @param filter     Apply a Butterworth low-pass filter
#' @param cutoff     Normalized cutoff frequency 
#' @param order      Order of the Butterworth filter (e.g. 2, 4, 6, 8).
#' @return A list with Parameters, Exp.Model (fit plot), RefLine.Model (residual plot), and Cor.Result (R2 / Spearman correlation table).
MonoExpModel <- function(data, variable, direction,
                         filter = TRUE, cutoff = 0.3, order = 2){
  
# Code for cleaning of data
# Removing rows with missing Time
  data <- data[!is.na(data$Time), ]
  
# Store original raw data for plotting later
  data$raw <- data[[variable]]
  data$.y  <- data[[variable]]
  
# Remove negative values
  neg_count <- sum(data$.y < 0, na.rm = TRUE)
  
  if(neg_count > 0){
    message(paste("Removed", neg_count, "negative values in", variable))
    data$.y[data$.y < 0] <- NA
  }
  
# Linear interpolation of NAs. NOTE: This is a linear interpolation. NOT Spline
  if(any(is.na(data$.y))){
    data$.y <- zoo::na.approx(
      data$.y,
      x = data$Time,
      na.rm = FALSE
    )
  }
  
# Fill any remaining NA values
  if(any(is.na(data$.y))){
    data$.y <- zoo::na.locf(data$.y, na.rm = FALSE)
    data$.y <- zoo::na.locf(data$.y, fromLast = TRUE)
  }
  
# Butterworth filter
  if(filter){
    bw <- signal::butter(order, cutoff, type = "low")
    data$.y <- signal::filtfilt(bw, data$.y)
  }
  
# Estimates 
  B_guess <- max(data$.y, na.rm = TRUE)
  
  threshold <- 0.05 * B_guess
  above_threshold <- data$Time[data$.y > threshold]
  
  TD_guess <- if(length(above_threshold) > 0){
    min(above_threshold, na.rm = TRUE)
  } else {
    0.05 * max(data$Time, na.rm = TRUE)
  }
  
  tau_guess <- 0.2 * (max(data$Time, na.rm = TRUE) - TD_guess)
  if(!is.finite(tau_guess) || tau_guess <= 0){
    tau_guess <- 0.2 * max(data$Time, na.rm = TRUE)
  }
  
  Start_vals <- list(
    B   = B_guess,
    tau = tau_guess,
    TD1 = TD_guess
  )
  
# Fit model 
  fit <- tryCatch({
    
    if(direction == 1){  # Rise
      minpack.lm::nlsLM(
        .y ~ B * (1 - exp(-pmax(Time - TD1, 0)/tau)),
        data = data,
        start = Start_vals,
        lower = c(B = 0, tau = 0, TD1 = 0),
        control = minpack.lm::nls.lm.control(maxiter = 300)
      )
    } else {  # Decay
      minpack.lm::nlsLM(
        .y ~ B * (1 - exp(-pmax(Time - TD1, 0)/tau)),
        data = data,
        start = Start_vals,
        control = minpack.lm::nls.lm.control(maxiter = 300)
      )
    }
    
  }, error = function(e){
    message("Model failed: ", e$message)
    return(NULL)
  })
  
  if(is.null(fit)){
    return(list(
      Parameters = NA,
      Exp.Model = ggplot() + ggtitle("Fit failed"),
      RefLine.Model = ggplot() + ggtitle("Fit failed"),
      Cor.Result = NA
    ))
  }
  
# Predictions
  data$Fit <- predict(fit)
  valid <- !is.na(data$raw) & !is.na(data$Fit)
  residuals <- data$raw[valid] - data$Fit[valid]
  SSres <- sum(residuals^2)
  SStot <- sum((data$raw[valid] - mean(data$raw[valid]))^2)
  
  if(SStot == 0){
    R2 <- NA
  } else {
    R2 <- 1 - SSres / SStot
  }
  
  cor_test <- cor.test(
    residuals,
    data$raw[valid],
    method = "spearman"
  )
  
# Plotting 
  model_plot <- ggplot(data, aes(x = Time)) +
    geom_point(aes(y = raw), color = "blue") +
    geom_line(aes(y = Fit), color = "red", linewidth = 1) +
    labs(y = variable) +
    theme_classic() +
    ggtitle(paste(variable, "Model Fit"))
  
  residual_plot <- ggplot(data[valid, ], aes(x = Time)) +
    geom_point(aes(y = residuals), color = "darkorange") +
    geom_hline(yintercept = mean(residuals), linetype = "dashed") +
    labs(y = "Residual") +
    theme_classic() +
    ggtitle(paste(variable, "Residuals"))
  
# Parameters & statistics
  # B: total amplitude
  # Tau: time constant (63% of B); 
  # TD: time delay;
  # MRT: mean response time (Tau + TD)
  params <- data.frame(
    Variable = variable,
    B   = coef(fit)[1],
    Tau = coef(fit)[2],
    TD  = coef(fit)[3],
    MRT = coef(fit)[2] + coef(fit)[3]
  )
  
  cor_results <- data.frame(
    R2  = R2,
    P   = cor_test$p.value,
    Rho = cor_test$estimate
  )
  
  list(
    Parameters = params,
    Exp.Model = model_plot,
    RefLine.Model = residual_plot,
    Cor.Result = cor_results
  )
}
