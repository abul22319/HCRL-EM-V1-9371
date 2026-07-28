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
                         n_comp = 1,
                         filter = TRUE, cutoff = 0.3, order = 2){

  # cleaning 
  data <- data[!is.na(data$Time), ]
  data$raw <- data[[variable]]
  data$.y  <- data[[variable]]

  neg_count <- sum(data$.y < 0, na.rm = TRUE)
  if(neg_count > 0){
    message(paste("Removed", neg_count, "negative values in", variable))
    data$.y[data$.y < 0] <- NA
  }
  if(any(is.na(data$.y))){
    data$.y <- zoo::na.approx(data$.y, x = data$Time, na.rm = FALSE)
  }
  if(any(is.na(data$.y))){
    data$.y <- zoo::na.locf(data$.y, na.rm = FALSE)
    data$.y <- zoo::na.locf(data$.y, fromLast = TRUE)
  }
  if(filter){
    bw <- signal::butter(order, cutoff, type = "low")
    data$.y <- signal::filtfilt(bw, data$.y)
  }

  n_comp <- as.numeric(n_comp)

  Tmax  <- max(data$Time, na.rm = TRUE)
  B_tot <- max(data$.y,   na.rm = TRUE)

  threshold <- 0.05 * B_tot
  above     <- data$Time[data$.y > threshold]
  TD_base   <- if(length(above) > 0) min(above, na.rm = TRUE) else 0.05 * Tmax

  tau_span <- Tmax - TD_base
  if(!is.finite(tau_span) || tau_span <= 0) tau_span <- Tmax

  # one B/tau/TD per component
  Start_vals <- list()
  for(i in seq_len(n_comp)){
    Start_vals[[paste0("B",   i)]] <- B_tot / n_comp
    Start_vals[[paste0("tau", i)]] <- 0.2 * tau_span * i / n_comp
    Start_vals[[paste0("TD",  i)]] <- TD_base
  }

  terms <- vapply(seq_len(n_comp), function(i){
    sprintf("B%d * (1 - exp(-pmax(Time - TD%d, 0)/tau%d))", i, i, i)
  }, character(1))
  model_formula <- as.formula(paste(".y ~", paste(terms, collapse = " + ")))

  # Fit
  fit <- tryCatch({
    if(direction == 1){
      lower_vec <- setNames(rep(0, length(Start_vals)), names(Start_vals))
      minpack.lm::nlsLM(
        model_formula, data = data, start = Start_vals,
        lower = lower_vec,
        control = minpack.lm::nls.lm.control(maxiter = 300)
      )
    } else {
      minpack.lm::nlsLM(
        model_formula, data = data, start = Start_vals,
        control = minpack.lm::nls.lm.control(maxiter = 300)
      )
    }
  }, error = function(e){
    message("Model failed: ", e$message)
    NULL
  })

  if(is.null(fit)){
    return(list(
      Parameters    = NA,
      Exp.Model     = ggplot() + ggtitle("Fit failed"),
      RefLine.Model = ggplot() + ggtitle("Fit failed"),
      Cor.Result    = NA
      )
    }
