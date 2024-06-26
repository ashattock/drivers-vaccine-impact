###########################################################
# COVERAGE
#
# All coverage related functionality in one place.
#
###########################################################

# ---------------------------------------------------------
# Parent function for preparing vaccine coverage data
# ---------------------------------------------------------
prepare_coverage = function() {
  
  # Extract coverage for VIMC pathogens
  vimc_dt = coverage_vimc()

  # However not every country is covered by VIMC for these pathogens
  vimc_countries_dt = vimc_dt %>%
    select(d_v_a_id, country, year, source) %>%
    arrange(d_v_a_id, country, year) %>%
    unique()

  # For other countries and years, extract coverage from WIISE database
  wiise_dt = coverage_wiise(vimc_countries_dt) %>%
    # Assume linear 1974-1980 scale up...
    linear_coverage_scaleup() %>%
    # Assume constant over most recent (post-COVID) years...
    constant_coverage_extapolation()
  
  # Incorporate non-routine SIA data (from WIISE)
  sia_dt = coverage_sia(vimc_countries_dt)  # See sia.R
  
  # Combine all coverage data sources
  source_dt = rbind(vimc_dt, wiise_dt, sia_dt) %>%
    # Deal with meningitis A special case...
    meningococcal_conjugate()
  
  # Sanity check that no zero entries remain
  if (any(source_dt$fvps <= 1e-10))
    stop("Trival coverage entries identified")
  
  # Summarise - assuming partially targeted SIAs - for all d-v-a
  everything_dt = source_dt %>%
    lazy_dt() %>%
    group_by(d_v_a_id, country, year, age) %>%
    summarise(fvps     = max(fvps),  # Essentially a placeholder until next calculation
              cohort   = mean(cohort), 
              coverage = 1 - prod(1 - coverage)) %>%  # Key assumption
    ungroup() %>%
    # Use combined coverage - unless FVPs already eclipses 100% (unlikely)
    mutate(fvps = pmax(cohort * coverage, fvps)) %>%
    as.data.table()
  
  # Subset of the d-v-a for which require EPI50 analytics
  coverage_dt = everything_dt %>%
    filter(d_v_a_id %in% table("d_v_a")$d_v_a_id)
  
  # Save this primary coverage datatable to file
  save_table(coverage_dt, "coverage")
  
  # Also save the supporting datatables to file
  save_table(source_dt,     "coverage_source")
  save_table(everything_dt, "coverage_everything")

  # ---- Data visualisation plots ----
  
  # Plot total number of FVP over time
  plot_total_fvps()
  
  # Coverage data density by age
  plot_coverage_age_density()
}

# ---------------------------------------------------------
# Extract coverage from VIMC outputs
# ---------------------------------------------------------
coverage_vimc = function() {
  
  message(" > Coverage data: VIMC")
  
  # Vaccines for which we'll use VIMC estimates
  d_v_a_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
    bind_rows(table("d_v_a_extern")) %>%
    select(d_v_a_id, disease, vaccine, activity)
  
  # Recode vaccine IDs for consistency
  vaccine_recode = c(
    hib3 = "hib", 
    pcv3 = "pcv", 
    rcv2 = "rubella")
  
  # Extract VIMC vaccine coverage data
  vimc_dt = fread(paste0(o$pth$input, "vimc_coverage.csv")) %>%
    select(disease, vaccine, activity = activity_type, country,
           gender, year, age, fvps_adjusted, cohort_size) %>%
    # Countries and timeframe of interest...
    filter(country %in% all_countries(), 
           year    %in% o$years) %>%
    # Recode disease and vaccines...
    mutate(disease = tolower(disease), 
           vaccine = tolower(vaccine), 
           vaccine = recode(vaccine, !!!vaccine_recode)) %>%
    # Rubella special case...
    mutate(activity = ifelse(
      test = disease == "rubella", 
      yes  = "all", 
      no   = activity)) %>%
    # Append d_v_a ID...
    inner_join(y  = d_v_a_dt, 
               by = c("disease", "vaccine", "activity")) %>%
    # Summarise where ...
    group_by(d_v_a_id, country, year, age) %>%
    summarise(fvps     = sum(fvps_adjusted),
              cohort   = sum(cohort_size),
              coverage = fvps / cohort) %>%
    ungroup() %>%
    # Tidy up...
    arrange(d_v_a_id, country, year, age) %>%
    mutate(source = "vimc") %>%
    as.data.table()
  
  return(vimc_dt)
}

# ---------------------------------------------------------
# Extract coverage from WIISE database
# ---------------------------------------------------------
coverage_wiise = function(vimc_countries_dt) {
  
  message(" > Coverage data: WIISE")
  
  # ---- Load data ----
  
  # File path for already-downloaded WIISE coverage data
  raw_file = paste0(o$pth$input, "wiise_coverage.csv")
  
  # If file has already been downloaded, read it now
  if (file.exists(raw_file)) {
    raw_dt = fread(raw_file)
    
  } else {  # Otherwise we'll need to download
    
    # Non-VIMC coverage taken from WIISE database
    raw_url = "https://whowiise.blob.core.windows.net/upload/coverage--2021.xlsx"
    raw_dt  = read_url_xls(raw_url, sheet = 1)
    
    # Save csv file locally for easy re-loading
    fwrite(raw_dt, file = raw_file)
  }
  
  # ---- Wrangle WIISE data ----
  
  # Routine activities (or 'all' for non-VIMC pathogens)
  d_v_a_dt = table("d_v_a") %>%
    filter(source != "extern") %>%
    bind_rows(table("d_v_a_extern")) %>%
    filter(activity %in% c("routine", "all")) %>%
    select(d_v_a_id, vaccine)
  
  # Parse 'interventions' into EPI50 vaccines
  reduced_dt = raw_dt %>%
    # Convert to lower case...
    setnames(names(.), tolower(names(.))) %>% 
    mutate_if(is.character, tolower) %>%
    # Reduce columns...
    select(intervention = antigen, country = code, 
           year, coverage, source = coverage_category)
  
  # Parse 'interventions' into EPI50 vaccines
  intervention_dt = reduced_dt %>%
    # Remove any unknown countries...
    mutate(country = toupper(country)) %>%
    filter(country %in% all_countries(), 
           year    %in% o$years) %>%
    # Convert coverage to proportion...
    mutate(coverage = coverage / 100) %>%
    filter(coverage > 0) %>%
    # Use WUENIC data as primary source...
    mutate(wuenic   = ifelse(source == "wuenic", coverage, NA), 
           coverage = ifelse(source != "wuenic", coverage, NA)) %>%
    # Compare against average of all other sources...
    lazy_dt() %>%
    group_by(country, intervention, year) %>%
    summarise(wuenic = mean(wuenic,   na.rm = TRUE),
              other  = mean(coverage, na.rm = TRUE)) %>%
    ungroup() %>%
    # Salvage coverage from non-WUENIC sources...
    mutate(wuenic = ifelse(is.nan(wuenic), other, wuenic)) %>%
    select(country, intervention, year, coverage = wuenic) %>%
    # Bound all non-trivial coverage values...
    mutate(coverage = pmin(coverage, o$max_coverage)) %>%
    filter(coverage > 0) %>%
    # Interpret 'intervention'...
    left_join(y  = table("vaccine_dict"), 
              by = "intervention", 
              relationship = "many-to-many") %>%
    filter(!is.na(vaccine)) %>%
    # Append d-v-a details...
    left_join(y  = d_v_a_dt, 
              by = "vaccine") %>%
    select(d_v_a_id, vaccine, intervention, country, year, coverage) %>%
    arrange(d_v_a_id, vaccine, intervention, country, year) %>%
    as.data.table()
  
  # Plot coverage value density by intervention ID
  # g = ggplot(intervention_dt) +
  #   aes(x = coverage,
  #       y = after_stat(count),
  #       colour = intervention,
  #       fill   = intervention) +
  #   geom_density(alpha = 0.2) +
  #   facet_wrap(~vaccine)
  
  # ---- Separately store global coverages ----
  
  # Global vaccine coverage according to WUENIC
  global_dt = reduced_dt %>%
    filter(country == "global", 
           source  == "wuenic") %>%
    select(-source) %>%
    # Interpret 'intervention'...
    left_join(y  = table("vaccine_dict"), 
              by = "intervention", 
              relationship = "many-to-many") %>%
    # Append d-v-a details...
    inner_join(y  = d_v_a_dt, 
               by = "vaccine") %>%
    select(d_v_a_id, vaccine, year, coverage) %>%
    arrange(d_v_a_id, vaccine, year) %>%
    # Coverage as a proportion...
    mutate(coverage = coverage / 100)
  
  # Save table in cache
  save_table(global_dt, "coverage_global")
  
  # ---- Calculate FVPs (non pregnancy vaccines) ----

  # Age at vaccination (deal with pregnancy vaccines after)
  age_dt = table("vaccine_age") %>%
    filter(age != "NA") %>%
    mutate(age = as.numeric(age))

  # Append age and calculate FVPs
  wiise_age_dt = intervention_dt %>%
    # Remove countries and years already covered by VIMC...
    left_join(y  = vimc_countries_dt, 
              by = c("d_v_a_id", "country", "year")) %>%
    filter(is.na(source)) %>%
    select(-intervention, -source) %>%
    # Append ages...
    left_join(y  = age_dt, 
              by = "vaccine") %>%
    filter(!is.na(age)) %>%
    # Calculate fully vaccinated people...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    mutate(sheduled_doses = coverage * pop) %>%
    calculate_fvps() %>%
    as.data.table()
  
  # ---- Calculate FVPs (pregnancy vaccines) ----
  
  # Age at vaccination for pregnancy vaccines
  age_birth_dt = table("vaccine_age") %>%
    left_join(y  = d_v_a_dt, 
              by = "vaccine") %>%
    filter(age == "NA") %>% 
    select(d_v_a_id) %>%
    expand_grid(table("wpp_fertility")) %>%
    # Remove trivial values...
    filter(fertility > 0) %>%
    group_by(country, year) %>%
    mutate(fertility = fertility / sum(fertility)) %>%
    ungroup() %>%
    as.data.table()
  
  # Append age and calculate FVPs
  wiise_pregnancy_dt = intervention_dt %>%
    select(-intervention) %>%
    # Reduce down to pregnancy vaccines...
    left_join(y  = age_dt, 
              by = "vaccine") %>%
    filter(is.na(age)) %>%
    # Coverage in this context is of newborns...
    mutate(age = 0) %>%
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    mutate(sheduled_doses = coverage * pop) %>%
    calculate_fvps() %>%
    # But we want FVPs in terms of mothers...
    select(-age, -cohort, -coverage) %>%
    left_join(y  = age_birth_dt,
              by = c("d_v_a_id", "country", "year"), 
              relationship = "many-to-many") %>%
    mutate(fvps = fvps * fertility) %>%
    # Append parental demographics...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    rename(cohort = pop) %>%
    # Calculate coverage (of mothers)...
    mutate(fvps = pmin(fvps, cohort * o$max_coverage),
           coverage = fvps / cohort) %>%
    # Tidy up...
    select(all_names(wiise_age_dt)) %>%
    as.data.table()
  
  # Combine into signle datatable
  wiise_dt = wiise_age_dt %>%
    rbind(wiise_pregnancy_dt) %>%
    arrange(d_v_a_id, country, year, age) %>%
    mutate(source = "wiise")
  
  return(wiise_dt)
}

# ---------------------------------------------------------
# FVP calculation considering boosters
# ---------------------------------------------------------
calculate_fvps = function(coverage_dt) {
  
  # NOTES: 
  #  - Using mean for pop as all values should all equal
  #  - Coverage bounded by o$max_coverage
  
  # For primary schedule, assume all new FVPs
  primary_dt = coverage_dt %>% 
    lazy_dt() %>%
    filter(!grepl("_bx$", vaccine)) %>%
    group_by(d_v_a_id, country, year, age) %>%
    summarise(fvps   = sum(sheduled_doses),  # Using sum
              cohort = mean(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # All booster doses
  booster_dt = coverage_dt %>% 
    filter(grepl("_bx$", vaccine))
  
  # Check whether trivial
  if (nrow(booster_dt) == 0) {
    booster_dt = NULL 
    
    } else {  # Otherwise, summarise...
    
    # For boosters, consecutive doses are for the same person
    booster_dt %<>%
      lazy_dt() %>%
      group_by(d_v_a_id, country, year, age) %>%
      summarise(fvps   = max(sheduled_doses),  # Using max
                cohort = mean(pop)) %>%
      ungroup() %>%
      as.data.table()
  }
    
  # Re-bind everything together and calculate coverage
  fvps_dt = rbind(primary_dt, booster_dt) %>%
    mutate(fvps = pmin(fvps, cohort * o$max_coverage),
           coverage = fvps / cohort)
  
  return(fvps_dt)
}

# ---------------------------------------------------------
# Assume a linear scale up prior to data start
# ---------------------------------------------------------
linear_coverage_scaleup = function(coverage_dt) {
  
  # Years we will scale up over
  scaleup_years = min(o$years) : (min(coverage_dt$year) - 1)
  
  # Income status in first year of data
  income_dt = coverage_dt %>%
    # Remove reference to FVPs, we'll recalculate...
    select(-fvps, -cohort) %>%
    # Append income status over time...
    left_join(y  = table("income_status"), 
              by = c("country", "year")) %>%
    # Non-trivial values from first year of data...
    filter(year == min(year))
  
  # Function to repeat trivialised coverage datatable for given year
  rep_fn = function(rep_year)
    income_dt %>% mutate(year = rep_year, coverage = NA)
  
  # For non-high-income countries, create blank scale up datatable
  template_dt = scaleup_years %>%
    # Repeat trivialised coverage datatable for each year
    lapply(rep_fn) %>%
    rbindlist() %>%
    rbind(income_dt) %>%
    arrange(d_v_a_id, country, age, year) %>%
    # Only interested in non-HIC...
    filter(income != "hic")
  
  # Set 1974 coverage to zero and linearly scale up to 1980
  if (o$pre_1980_assumption == "linear") {
    
    # NOTE: A conservative assumption
    scaleup_dt = template_dt %>%
      # Start at zero coverage...
      mutate(coverage = ifelse(
        test = year == min(scaleup_years), 
        yes  = 0, 
        no   = coverage)) %>%
      # Linearly interpolate from zero to 1980 coverage...
      group_by(d_v_a_id, country, age) %>%
      mutate(coverage = na_interpolation(coverage)) %>%
      ungroup() %>%
      as.data.table()
  }
  
  # Alternatively assume constant over this period
  if (o$pre_1980_assumption == "constant") {
    
    # NOTE: An ambitious assumption
    scaleup_dt = template_dt %>%
      group_by(d_v_a_id, country, age) %>%
      fill(coverage, .direction = "up") %>%
      ungroup() %>%
      as.data.table()
  }
  
  # Append cohort size and calculate FVPs
  scaleup_dt %<>%
    # Remove 1980 value to avoid repetition...
    filter(year %in% scaleup_years, 
           coverage > 0) %>%
    # Calculate FVPs...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    rename(cohort = pop) %>%
    mutate(fvps = cohort * coverage) %>%
    # Tidy up...
    select(all_names(coverage_dt)) %>%
    arrange(d_v_a_id, country, year, age) %>%
    as.data.table()
  
  # For high-income countries, assume constant over this period
  constant_dt = income_dt %>%
    filter(income == "hic") %>%
    # KEY ASSUMPTION: Repeat coverage for early years...
    select(-year) %>%
    expand_grid(year = scaleup_years) %>%
    # Append cohort size and calculate FVPs...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    rename(cohort = pop) %>%
    mutate(fvps = cohort * coverage) %>%
    # Tidy up...
    select(all_names(coverage_dt)) %>%
    arrange(d_v_a_id, country, year, age) %>%
    as.data.table()
  
  # Bind these two datatables into coverage
  coverage_dt %<>%
    rbind(scaleup_dt) %>%
    rbind(constant_dt) %>%
    arrange(d_v_a_id, country, year, age)
  
  return(coverage_dt)
}

# ---------------------------------------------------------
# Assume constant coverage over most recent years
# ---------------------------------------------------------
constant_coverage_extapolation = function(coverage_dt) {
  
  # Extrapolate coverage data from most recent year
  extrap_dt = coverage_dt %>%
    # Remove reference to FVPs, we'll recalculate...
    select(-fvps, -cohort) %>%
    # Years from which to extrapolate (3 years with the past 5)...
    filter(year >= max(o$years) - 5) %>%
    group_by(d_v_a_id, country, age) %>%
    slice_max(year, n = 3, with_ties = FALSE) %>%
    ungroup() %>%
    # Mean coverage over these recent years...
    group_by(d_v_a_id, country, age) %>%
    summarise(year     = max(year) + 1, 
              coverage = mean(coverage)) %>%
    ungroup() %>%
    # KEY ASSUMPTION: Repeat coverage for most recent years...
    expand_grid(extrap_year = o$years) %>%
    group_by(d_v_a_id, country) %>%
    filter(extrap_year >= year) %>%
    ungroup() %>%
    select(d_v_a_id, country, age,
           year = extrap_year, coverage) %>%
    # Append cohort size and calculate FVPs...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    rename(cohort = pop) %>%
    mutate(fvps   = cohort * coverage, 
           source = "wiise") %>%
    # Tidy up...
    select(all_names(coverage_dt)) %>%
    arrange(d_v_a_id, country, year, age) %>%
    as.data.table()
  
  # Bind these two datatables into coverage
  coverage_dt %<>%
    rbind(extrap_dt) %>%
    arrange(d_v_a_id, country, year, age)
  
  return(coverage_dt)
}

# ---------------------------------------------------------
# Remove effects on men_conj vaccine in locations without menA burden
# ---------------------------------------------------------
meningococcal_conjugate = function(coverage_dt) {
  
  # Meningitis A d-v-a IDs
  mena_id = table("d_v_a") %>%
    filter(disease == "mena") %>%
    pull(d_v_a_id)
  
  # Meningitis A belt according to VIMC
  mena_countries = coverage_dt %>%
    filter(d_v_a_id %in% mena_id, 
           source == "vimc") %>%
    pull(country) %>%
    unique()
  
  # Remove coverage for countries outside of the Men A belt
  mena_belt_dt = coverage_dt %>%
    filter(d_v_a_id %in% mena_id, 
           country  %in% mena_countries)
  
  # Update coverage datatable
  updated_dt = coverage_dt %>%
    filter(!d_v_a_id %in% mena_id) %>%
    rbind(mena_belt_dt) %>%
    arrange(d_v_a_id, country, year, age)
  
  return(updated_dt)
}

