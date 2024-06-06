###########################################################
# PREPARE
#
# Prepare various data sources for use throughout the analysis.
# The idea is that this process needs to be done only once, 
# and that the prepared inputs are streamlined for quick loading.
#
###########################################################

# ---------------------------------------------------------
# Parent function for preparing model inputs from various data sources
# ---------------------------------------------------------
run_prepare = function() {
  
  # Only continue if specified by run_module
  if (!is.element(1, o$run_module)) return()
  
  message("* Preparing input data")
  
  # Convert config yaml files to datatables
  prepare_config_tables()
  
  # Streamline VIMC impact estimates for quick loading
  prepare_vimc_estimates()
  
  # Prepare country income status classification over time
  prepare_income_status()
  
  # Prepare demography-related estimates from WPP
  prepare_demography()
  
  # Prepare all covariates for regression modelling
  prepare_covariates()  # See covariates.R
  
  # Prepare historical vaccine coverage
  prepare_coverage()  # See coverage.R
}

# ---------------------------------------------------------
# Convert config yaml files to datatables 
# ---------------------------------------------------------
prepare_config_tables = function() {
  
  # NOTE: Convert from yaml (/config) to rds (/tables) for fast loading
  
  # List of config yaml files to convert
  config_files = o$pth$config %>%
    list.files(pattern = ".+\\.yaml$") %>%
    str_remove(".yaml$") %>%
    setdiff("general")
  
  # Iterate through these files
  for (file in config_files) {
    
    # Load the yaml file
    yaml_file = paste0("config/", file, ".yaml")
    yaml_data = read_yaml(yaml_file)$table
    
    # Convert to datatable
    config_dt = yaml_data %>%
      lapply(as.data.table) %>%
      rbindlist(fill = TRUE)
    
    # Save in tables cache
    save_table(config_dt, file)
  }
}

# ---------------------------------------------------------
# Streamline VIMC impact estimates for quick loading
# ---------------------------------------------------------
prepare_vimc_estimates = function() {
  
  message(" > VIMC estimates")
  
  # All diseases to load VIMC outcomes for
  vimc_info = table("d_v_a") %>%
    filter(source == "vimc") %>%
    select(disease) %>%
    unique() %>%
    left_join(y  = table("disease_name"), 
              by = "disease")
  
  # Initiate list to store outcomes
  vimc_list = list()
  
  # Iterate through diseases
  for (i in seq_row(vimc_info)) {
    
    # Disease ID and associated full name
    id   = vimc_info[i]$disease
    name = vimc_info[i]$disease_name
    
    message("  - ", name)
    
    # Load VIMC impact estimates for this disease
    vimc_list[[id]] = read_rds("vimc", id) %>%
      lazy_dt() %>%
      pivot_longer(cols = ends_with("impact"),
                   names_to = "vaccine") %>%
      replace_na(list(value = 0)) %>%
      # Intrepet disease, vaccine, and activity...
      mutate(disease  = tolower(disease),
             vaccine  = str_remove(vaccine, "_impact"),
             activity = ifelse(
               test = vaccine %in% c("routine", "campaign"),
               yes  = vaccine,
               no   = "routine")) %>%
      # Tidy up...
      select(disease, vaccine, activity, country,
             year, age, metric = outcome, value) %>%
      as.data.table()
  }
  
  # Squash results into single datatable
  vimc_dt = rbindlist(vimc_list) %>%
    # Interpret activity...
    mutate(vaccine = ifelse(
      test = vaccine %in% c("routine", "campaign"),
      yes  = disease,
      no   = vaccine)) %>%
    # Deal with rubella special case...
    mutate(is_all   = vaccine == "combined",
           vaccine  = ifelse(is_all, disease, vaccine),
           activity = ifelse(is_all, "all", activity)) %>%
    # Wide format of metrics...
    mutate(metric = paste1(metric, "averted")) %>%
    pivot_wider(names_from = metric) %>%
    # Append d-v-a ID...
    left_join(y  = table("d_v_a"),
              by = c("disease", "vaccine", "activity")) %>%
    # Tidy up...
    select(d_v_a_id, country, year, age,
           deaths_averted, dalys_averted) %>%
    arrange(d_v_a_id, country, year, age) %>%
    as.data.table()
  
  # Save in tables cache
  save_table(vimc_dt, "vimc_estimates")
}

# ---------------------------------------------------------
# Prepare country income status classification over time
# ---------------------------------------------------------
prepare_income_status = function() {
  
  message(" > Income status")
  
  # Path to data file
  #
  # SOURCE: https://datacatalogfiles.worldbank.org/ddh-published/0037712/
  #         DR0090755/CLASS.xlsx?versionId=2023-11-16T18:35:30.5758473Z
  #
  # Alternatively, download 'Historical classification by income' Excel file from: 
  # datacatalog.worldbank.org/search/dataset/0037712/World-Development-Indicators
  file = paste0(o$pth$input, "worldbank_income_status.csv")
  
  # Full country-year combination
  full_dt = expand_grid(
    country = all_countries(), 
    year    = o$years) %>%
    as.data.table()
  
  # Load and format country income status over time
  income_dt = fread(file, header = TRUE) %>%
    # Countries of interest...
    filter(country %in% all_countries()) %>%
    select(-country_name) %>%
    # Convert to tidy format...
    pivot_longer(cols = -country, 
                 names_to  = "year", 
                 values_to = "income") %>%
    # Country with all full country-year combo...
    mutate(year = as.integer(year)) %>%
    full_join(y  = full_dt, 
              by = c("country", "year")) %>%
    arrange(country, year) %>%
    # Fill missing data with pro/preceding value...
    mutate(income = ifelse(income == "", NA, income)) %>%
    group_by(country) %>%
    fill(income, .direction = "downup") %>%
    ungroup() %>%
    # Niue and Cook Islands missing, both are HIC...
    replace_na(list(income = "H")) %>%
    mutate(income = paste0(tolower(income), "ic")) %>%
    as.data.table()
  
  # Save in tables cache
  save_table(income_dt, "income_status")
}

# ---------------------------------------------------------
# Prepare demography-related estimates from WPP
# ---------------------------------------------------------
prepare_demography = function() {
  
  message(" > Demography data")
  
  # Function to apply element-wise scaler to data
  scaler_fn = function(m) {
    
    # Population scaling: by country, year, and age
    if (grepl("pop", m$scale))
      scaler_dt = setnames(
        x   = table("wpp_pop"), 
        old = "pop", 
        new = "scaler")
    
    # Numeric values: simple repitition
    if (grepl("^[0-9,\\.]+$", m$scale))
      scaler_dt = expand_grid(
        country = all_countries(), 
        year    = o$years,
        age     = if (m$age) o$ages else NA, 
        scaler  = as.numeric(m$scale)) %>%
        as.data.table()
    
    return(scaler_dt)
  }
  
  # Details of WPP metrics to load
  wpp_metrics = table("wpp_dict") %>%
    mutate(metric = fct_inorder(metric)) %>%
    split(.$metric)
  
  # Iterate through metrics to load
  for (metric in names(wpp_metrics)) {
    m = as.list(wpp_metrics[[metric]])
    
    message("  - ", metric)
    
    # Past and future in separate data sets
    if (m$proj == TRUE) {
      
      # Age-disaggregation specified in data set file name
      age = ifelse(m$age, "Age", "")
      
      # Names of WPP2022 data files to load
      past   = paste0(m$file,         age, o$pop_bin, "dt")
      future = paste0(m$file, "proj", age, o$pop_bin, "dt") 
      
      # Load pop data from WPP github package
      data_list = data_package(past, future, package = "wpp2022")
    }
    
    # Past and future combined into single data set
    if (m$proj == FALSE) {
      
      # Name of WPP2022 data file - history and projection in one
      all_time = paste0(m$file, o$pop_bin, "dt")
      
      # Load data from WPP github package
      data_list = data_package(all_time, package = "wpp2022")
    }
    
    # Combine (extended) past and future data
    data_dt = rbindlist(data_list, fill = TRUE) %>%
      {if (!m$age) mutate(., age = NA) else .} %>%
      # Select countries of interest...
      inner_join(y  = table("country"),  
                 by = "country_code") %>%
      select(country, year, age, value = !!m$var) %>%
      # Shift year by one (see github.com/PPgp/wpp2022 for details)...
      mutate(year = as.integer(year) + 1) %>%
      filter(year %in% o$years) %>%
      # Scale metrics...
      left_join(y = scaler_fn(m), 
                by = c("country", "year", "age")) %>%
      mutate(value = value * scaler) %>%
      select(-scaler) %>%
      # Tidy up...
      rename(!!metric := value) %>%
      arrange(country, year, age)
    
    # Save in tables cache
    save_table(data_dt, paste1("wpp", metric))
  }
}

# ---------------------------------------------------------
# Simple wrapper to load all countries
# ---------------------------------------------------------
all_countries = function(as_dt = FALSE) {
  
  # Pull all countries defined in config file
  countries = table("country")$country
  
  # Convert to simple datatable if desired
  if (as_dt == TRUE)
    countries = data.table(country = countries)
  
  return(countries)
}

# ---------------------------------------------------------
# Simple wrapper to load all regions
# ---------------------------------------------------------
all_regions = function() {
  
  # Pull all regions defined in config file
  regions = table("region_dict")$region
  
  return(regions)
}

# ---------------------------------------------------------
# Save table in cache directory for quick loading
# ---------------------------------------------------------
save_table = function(x, table) {
  
  # Save table in tables cache directory
  save_rds(x, "tables", table, "table")
}

# ---------------------------------------------------------
# Load and return cached datatable
# ---------------------------------------------------------
table = function(table) {
  
  # Construct file path
  file = paste0(o$pth$tables, table, "_table.rds")
  
  # Throw an error if this file doesn't exist
  if (!file.exists(file))
    stop("Table ", table, " has not been cached")
  
  # Load rds file
  y = read_rds(file)
  
  return(y)
}

