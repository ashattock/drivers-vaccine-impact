###########################################################
# UNCERTAINTY
#
# xxxxxxx
#
###########################################################

# ---------------------------------------------------------
# Parent function for generating uncertainty
# Called by: launch.R, main.R (and other launch-style scripts)
# ---------------------------------------------------------
run_uncertainty = function() {
  
  # Only continue if specified by do_step
  if (!is.element(3, o$do_step)) return()
  
  message("* Generating uncertainty")
  
  # Load impact factors calculated in step 2
  impact_dt       = try_load(o$pth$impact_factors, "impact_dt")
  scenario_impact = try_load(o$pth$impact_factors, "scenario_impact")
  
  # ---- Draws for VIMC diseases ----
  
  message(" - VIMC diseases")
  
  # TODO: Look into raking, and tidy up this section...
  
  # NOTE: Do we not need to consider activity_type here? As raking is vaccine-activity combo.
  raking_dt = summarize_raking(impact_dt) %>%
    select(vaccine, cv)
  
  # Append or calculate (using cv) standard deviation for all VIMC diseases
  #
  # NOTE: vimc_ui generated by data-raw/prep_vimc_ui.R
  vimc_dt = vimc_ui %>%
    right_join(y  = scenario_impact,  # Join impact factor
               by = c("disease", "country", "year")) %>%
    filter(disease %in% o$disease$vimc) %>%  # VIMC diseases only
    left_join(raking_dt, by = "vaccine") %>%
    mutate(sd = ifelse(is.na(sd), deaths_averted * cv, sd)) %>%  # ?? Imputing missing values (for non-VIMC countries) ??
    select(disease, vaccine, activity_type, country,
           impact_factor, year, age, fvps, deaths_averted, sd)
  
  # Generate uncertainty for VIMC diseases
  vimc_draws = generate_vimc_uncertainty(vimc_dt)
  
  # ---- Draws for non-VIMC diseases ----
  
  message(" - GBD diseases")
  
  # Generate draws for non-VIMC diseases using CI of initial vaccine efficacy
  #
  # NOTE: gbd_efficacy generated by data-raw/gbd_efficacy.R
  gbd_dt = gbd_efficacy %>%
    filter(disease %in% o$disease$gbd) %>%
    select(-source)
  
  # Generate uncertainty for GBD diseases
  gbd_draws = generate_gbd_uncertainty(gbd_dt, scenario_impact)
  
  # ---- Combine and align to modelled mean ----
  
  message(" - Realigning means")
  
  # Bind and shift all draws such that we don't move the modelled mean
  draws_dt = rbind(gbd_draws, vimc_draws) %>%                  # Concatenate all diseases
    mutate(draw_mean = rowMeans(across(starts_with("draw"))),  # Mean of all draws per row
           draw_diff = deaths_averted - draw_mean,             # Difference to modelled mean
           across(starts_with("draw"), ~. + draw_diff)) %>%    # Shift all draws to align means
    select(-draw_mean, -draw_diff)
  
  # TODO: Throw error if any NAs or <= 0 ??
  
  # Save datatable to file
  save_file(draws_dt, o$pth$uncertainty, "draws")
  
  # ---- Diagnostic figures ----
  
  # Check flag
  if (o$plot_diagnostics) {
    
    message(" - Plotting diagnostics")
    
    # Plot parameters of fitted beta distribution to vaccine efficacy (GBD diseases only)
    fig_name = "Uncertainty distributions - GBD diseases"
    plot_gbd_uncertainty_dist(fig_name)  # See plotting.R
    
    # Plot parameters of fitted beta distribution to vaccine efficacy (GBD diseases only)
    fig_name = "Uncertainty distribution fit - GBD diseases"
    plot_gbd_uncertainty_fit(fig_name)  # See plotting.R
    
    # Plot uncertainty draws for all diseases
    fig_name = "Uncertainty draws - All diseases"
    plot_draws(fig_name)  # See plotting.R
    
    # Plot annual totals to check alignment of means
    fig_name = "Uncertainty bounds - Annual total"
    plot_annual_total(fig_name)  # See plotting.R
  }
}

# ---------------------------------------------------------
# Sample VIMC uncertainty using mean and sd from provided results
# Called by: run_uncertainty()
# ---------------------------------------------------------
generate_vimc_uncertainty = function (vimc_dt) {
  
  # Sample draws from normal distribution
  draws_dt = rnorm %>%  # Sample using normal distribution
    mapply(n    = o$n_draws, 
           mean = vimc_dt$deaths_averted, 
           sd   = vimc_dt$sd) %>%
    t() %>%  # One draw per column
    as_named_dt(paste0("draw_", 1 : o$n_draws))
  
  # Append draws to details 
  draws_wide = cbind(vimc_dt, draws_dt) %>%
    select(-sd)
  
  return(draws_wide)
}

# ---------------------------------------------------------
# Sample GBD uncertainty using initial efficacy confidence intervals
# Called by: run_uncertainty()
# ---------------------------------------------------------
generate_gbd_uncertainty = function(gbd_dt, scenario_impact) {
  
  # Initiate optimal matrix (diseases x beta distribution parameters)
  opt_pars = matrix(NA, nrow = nrow(gbd_dt), ncol = 2)
  
  # ---- Optimise parameters for each disease ----
  
  # Loop through diseases
  for (i in 1 : nrow(gbd_dt)) {
    
    # Vaccine efficacy details
    v = efficacy[i, .(mean, lower, upper)]
    
    # Determine optimal beta parameters to fit vaccine efficacy
    opt_result = optim(par    = c(1, 1),     # Starting point
                       fn     = gbd_obj_fn,  # Objective function
                       data   = as.list(v),  # Additonal arguments for obj_fn
                       lower  = o$par_lower, 
                       upper  = o$par_upper,
                       method = "L-BFGS-B")
    
    # Store best fitting parameters
    opt_pars[i, ] = exp(opt_result$par)
  }
  
  # Convert optimal results to datatable
  beta_pars = opt_pars %>%
    as_named_dt(c("p1", "p2")) %>%
    mutate(disease = gbd_dt$disease, 
           vaccine = gbd_dt$vaccine, 
           .before = 1)
  
  # Save file - primarily for diagnostic figures
  save_file(beta_pars, o$pth$uncertainty, "gbd_beta_pars")
  
  # ---- Draw samples from optimal beta distribution ----
  
  # Draw samples from beta distribution using these optimal parameters
  draws_dt = beta_pars %>%
    # Sample from beta distribution...
    split(rownames(.)) %>%
    lapply(gbd_draw_fn, n = o$n_draws) %>%
    as_named_dt(efficacy$disease) %>%
    # Melt to long format...
    mutate(draw = paste0("draw_", 1 : o$n_draws)) %>%
    pivot_longer(cols = -draw, 
                 names_to  = "disease", 
                 values_to = "scaler") %>%
    select(disease, draw, scaler) %>%
    arrange(disease) %>%
    as.data.table()
  
  # Join with reference results and apply scaler
  draws_wide = draws_dt %>%
    left_join(scenario_impact, by = "disease") %>%
    mutate(deaths_averted_draw = deaths_averted * scaler) %>%
    select(-scaler) %>%
    pivot_wider(names_from  = draw, 
                values_from = deaths_averted_draw) %>%
    as.data.table()
  
  return(draws_wide)
}

# ---------------------------------------------------------
# Objective function to minimise
# ---------------------------------------------------------
gbd_obj_fn = function(par, data) {
  
  # Extract fitting parameters
  a = exp(par[1])
  b = exp(par[2])
  
  # Error between means
  mean_diff = abs(a / (a + b) - data$mean)
  
  # Error between bounds
  l_diff = abs(qbeta(0.05, a, b) - data$lower)
  u_diff = abs(qbeta(0.95, a, b) - data$upper)
  
  # Add weight to any error from the means
  y = sum(c(20 * mean_diff, l_diff, u_diff))
  
  return(y)
}

# ---------------------------------------------------------
# Draw samples for beta distribution (centered to mean)
# ---------------------------------------------------------
gbd_draw_fn = function(x, n) {
  
  # NOTE: alpha / (alpha + beta) is the mean of a beta distribution
  rbeta(n, x$p1, x$p2) / (x$p1 / (x$p1 + x$p2))
}
