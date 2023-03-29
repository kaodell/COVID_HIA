# COVID_HIA_README
This repository contains python3 scripts written by Katelyn O’Dell for calculating the health impact assessment in He et al., "COVID-19 perturbation on US air quality and health impact assessment" in prep for PNAS. If you have any questions about these scripts, want to use them, or find any errors please contact Katelyn O’Dell at the email address listed on her github home page.

Last updated: 03.29.2023

Contents:
- prep_pop_data_Jian.py, creates regridded population data for use in the Health Impact Assessment (HIA)
- grid_baseline_mort_Jian.py, takes county-level baseline mortality data and grids them to the WRF-Chem grid
- create_daily_files_Jian.py, takes monthly files from Jian WRF-Chem output and converts to daily files to use in the local daily average code
- create_loc_daily_average.py, takes the output of create_daily_files_Jian and calculates daily means (or mda8 for ozone) of pollutants using local time 
- run_HIA_gridded.py, calculates attributable mortality for annual average PM2.5 and O3 in the BAU and COV cases using a health impact assessment
- mk_figures_4Jian.py, creates figure 4 in the manuscript
- ODell_udf_Jian.py, contains functions used in other scripts 
- covid_hia_pyenv.txt, python3 environment used to run these scripts

## To recreate HIA analysis and figures in “insert Jian title” run the scripts in this folder in the following order:

### Step 1: Prepare datasets

#### 1a: run prep_pop_data_Jian.py 
to create regridded population data for use in the Health Impact Assessment (HIA)\
Inputs:
- NASA SECDAC GPW population density v4.11 at 2.5 arcminutes
- USCB 2018 national US shapefile at 5m
- Grid from WRF-Chem output from Dr. Jian He

Outputs:
- netCDF file containing regridded population

#### 1b: run grid_baseline_mort_Jian.py 
to take county-level baseline mortality data and grid them to the WRF-Chem grid. This script takes about 15 minutes to run on my local machine for the 12 km grid.
Inputs:
- Grid from WRF-Chem output from Dr. Jian He
- CDC WONDER county and state-level all cause baseline mortality rates for all ages, mean for 2015-2019. Details on selections made to download this data and access dates are below.
- 2019 TIGER/Line county shape files
Outputs:
- A csv containing flattened grid cell lats, lons and assigned baseline mortality rates in units of deaths per 100,000

#### 1c: run create_daily_files_Jian.py 
to take monthly files from Jian WRF-Chem output and converts to daily files to use in the local daily average code. This script takes about 30 minutes to run on my local machine
Inputs:
- Monthly  WRF-Chem output from Jian He
Outputs:
- Daily WRF-Chem files

#### 1d: run create_loc_daily_average.py
Takes the output of create_daily_files_Jian and calculates daily means (or mda8 for ozone) of pollutants using local time. For pm2.5, the code calculates wet pm2.5 mass at 35% relative humidity. This code has to be re-run for the desired pollutants (o3 and pm2.5) as well as the two simulation cases (bau and cov).
Inputs: 
- Daily WRF-Chem files output from create_daily_files_Jian.py
Outputs:
- Single file (for each run) containing local daily averages for each grid cell during each day of the specified time period.

### Step 2: run health impact assessment
#### run_HIA_gridded.py
Calculated attributable mortality for annual average PM2.5 and O3 in the BAU and COV cases using a health impact assessment. This code has to be run separately for PM2.5 an dO3.
Inputs:
- Annual files of local daily mean PM2.5 and local mda8 O3 output from create_loc_daily_average.py
- Regridded population output from prep_pop_data_Jian.py
- Gridded baseline mortality rates from grid_baseline_mort_Jian.py
Outputs
- npz file of annual mean pollutant concentrations and attributable mortality for both simulation cases

### Step 3: create figures
#### mk_figures_4Jian.py
Create figure 4 in the manuscript.
Inputs:
- O3 and pm2.5 hia results from run_HIA_gridded.py
- Regridded population from prep_pop_data_Jian.py
- USCB national US shapefile at 5m
Outputs:
- Figure 4 in the manuscript

### Additional codes used/information needed 
#### ODell_udf_Jian.py
Contains functions used in other scripts 

#### covid_hia_pyenv.txt
Python3 environment used to run these scripts.

### Details on dataset access dates and download parameters

CDC WONDER baseline mortality rates, state level
"Dataset: Underlying Cause of Death, 1999-2020"
"Query Parameters:"
"Year/Month: 2015; 2016; 2017; 2018; 2019"
"Group By: State"
"Show Totals: True"
"Show Zero Values: True"
"Show Suppressed: True"
"Calculate Rates Per: 100,000"
"Rate Options: Default intercensal populations for years 2001-2009 (except Infant Age Groups)"
"---"
"Help: See http://wonder.cdc.gov/wonder/help/ucd.html for more information."
"---"
"Query Date: Sep 20, 2022 9:58:06 AM"

CDC WONDER baseline mortality rates, county level
"Dataset: Underlying Cause of Death, 1999-2020"
"Query Parameters:"
"Year/Month: 2015; 2016; 2017; 2018; 2019"
"Group By: County"
"Show Totals: True"
"Show Zero Values: True"
"Show Suppressed: True"
"Calculate Rates Per: 100,000"
"Rate Options: Default intercensal populations for years 2001-2009 (except Infant Age Groups)"
"---"
"Help: See http://wonder.cdc.gov/wonder/help/ucd.html for more information."
"---"
"Query Date: Sep 20, 2022 9:16:35 AM"
