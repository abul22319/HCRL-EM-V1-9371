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

  Tmax <- max(data$Time, na.rm = TRUE)


  n_pts    <- length(data$.y)
  k        <- max(3, floor(0.10 * n_pts))          # points averaged at each end
  baseline <- mean(head(data$.y, k), na.rm = TRUE) # starting level
  plateau  <- mean(tail(data$.y, k), na.rm = TRUE) # ending level
  amp      <- plateau - baseline                   # SIGNED total amplitude
  if(!is.finite(amp) || amp == 0) amp <- diff(range(data$.y, na.rm = TRUE))


  thr     <- baseline + 0.05 * amp
  crossed <- if(amp >= 0) which(data$.y > thr) else which(data$.y < thr)
  TD_base <- if(length(crossed) > 0) data$Time[min(crossed)] else 0.05 * Tmax

  tau_span <- Tmax - TD_base
  if(!is.finite(tau_span) || tau_span <= 0) tau_span <- Tmax



  Start_vals <- list()
  for(i in seq_len(n_comp)){
    Start_vals[[paste0("B",   i)]] <- amp / n_comp
    Start_vals[[paste0("tau", i)]] <- 0.2 * tau_span * i / n_comp
    Start_vals[[paste0("TD",  i)]] <- TD_base + (i - 1) * tau_span / n_comp
  }

  terms <- vapply(seq_len(n_comp), function(i){
    sprintf("B%d * (1 - exp(-pmax(Time - TD%d, 0)/tau%d))", i, i, i)
  }, character(1))
  model_formula <- as.formula(paste(".y ~", paste(terms, collapse = " + ")))


  lower_vec <- setNames(numeric(length(Start_vals)),      names(Start_vals))
  upper_vec <- setNames(rep(Inf, length(Start_vals)),     names(Start_vals))
  for(i in seq_len(n_comp)){
    lower_vec[[paste0("tau", i)]] <- 1e-6
    lower_vec[[paste0("TD",  i)]] <- 0
    upper_vec[[paste0("TD",  i)]] <- Tmax
    lower_vec[[paste0("B",   i)]] <- if(direction == 1) 0 else -Inf
  }

  fit <- tryCatch({
    minpack.lm::nlsLM(
      model_formula, data = data, start = Start_vals,
      lower = lower_vec, upper = upper_vec,
      control = minpack.lm::nls.lm.control(maxiter = 300)
    )
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
    ))
  }

  # Predictions / R2 / correlation
  data$Fit  <- predict(fit)
  valid     <- !is.na(data$raw) & !is.na(data$Fit)
  residuals <- data$raw[valid] - data$Fit[valid]
  SSres <- sum(residuals^2)
  SStot <- sum((data$raw[valid] - mean(data$raw[valid]))^2)
  R2    <- if(SStot == 0) NA else 1 - SSres / SStot

  cor_test <- cor.test(residuals, data$raw[valid], method = "spearman")

  # Plots
  model_plot <- ggplot(data, aes(x = Time)) +
    geom_point(aes(y = raw), color = "blue") +
    geom_line(aes(y = Fit), color = "red", linewidth = 1) +
    labs(y = variable) + theme_classic() +
    ggtitle(paste(variable, "Model Fit —", n_comp, "component(s)"))

  residual_plot <- ggplot(data[valid, ], aes(x = Time)) +
    geom_point(aes(y = residuals), color = "darkorange") +
    geom_hline(yintercept = mean(residuals), linetype = "dashed") +
    labs(y = "Residual") + theme_classic() +
    ggtitle(paste(variable, "Residuals"))

  # Parameter table adapts to n_comp
  cf     <- coef(fit)
  params <- data.frame(Variable = variable, Components = n_comp)
  for(i in seq_len(n_comp)){
    params[[paste0("B",   i)]] <- cf[[paste0("B",   i)]]
    params[[paste0("Tau", i)]] <- cf[[paste0("tau", i)]]
    params[[paste0("TD",  i)]] <- cf[[paste0("TD",  i)]]
    params[[paste0("MRT", i)]] <- cf[[paste0("tau", i)]] + cf[[paste0("TD", i)]]
  }

  cor_results <- data.frame(
    R2  = R2,
    P   = cor_test$p.value,
    Rho = cor_test$estimate
  )

  list(
    Parameters    = params,
    Exp.Model     = model_plot,
    RefLine.Model = residual_plot,
    Cor.Result    = cor_results
  )
}
