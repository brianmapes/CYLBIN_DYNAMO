CYLBIN_output_SMARTR_legs_all.sav / CYLBIN_output_SMARTR_legs_all.nc


This file contains the output from the CYLBIN processing for the SMARTR data for the DYNAMO period

The descriptions of the variable stored in the output file is given below:


timestrings: time string in format YYYYDDMMHH

range: radius from the center where the divergence estimates are made, in km [size: 12 ranges]

Pressure: pressure in mb [size: 19 levels]

Divergence: divergence estimates in per sec [size: #hours x 12 x 19]

Divergence_stdev: uncertainty in the divergence estimates in per sec units [size: #hours x 12 x 19]
 
U: zonal winds estimated at each range radius, in m/s [size: #hours x 19]

V: meridional winds estimated at each range radius, in m/s [size: #hours x 19]

U_stdev: uncertainty in zonal winds zonal winds estimated at each range radius, in m/s [size: #hours x 19]

V_stdev: uncertainty in meridional winds estimated at each range radius, in m/s [size: #hours x 19]


chi_uv: chi-square coefficient of VAD fitted winds to the data at each range radius, in m/s [size: #hours x 19]

chi_div: chi-square coefficients of VAD divergence estimates in per sec [size: #hours x 12 x 19]
 
 
 EchoCover_0db: Fractional coverage of reflectivity bins containing values > 0 dBZ at a given range radius [# hours x 12 x 19]
 
 
 EchoCover_15db: Fractional coverage of reflectivity bins containing values > 15 dBZ at a given range radius [# hours x 12 x 19]
 
 EchoCover_30db: Fractional coverage of reflectivity bins containing values > 30 dBZ at a given range radius [# hours x 12 x 19]
 
 Arearain_GATEZR: area rain averaged with in the given range radius using GATE Z-R relation[# hours x 12] 

 Arearain_localZR; area rain averaged with in the given range radius using GATE Z-R relationship[# hours x 12] 

 stratiform_fraction: stratiform rain fraction with in the given range radius [# hours x 12] 
 

 All_data_matrix /v01: Master array containing 29 variables (listed below) stored in detail in a  single array: [size: 29 x 12 x 19 x # hours]

All_data_matrix = ['cover total'          ,'dbz (db)'            ,'Ze (mm^6/m^3)'       ,$
    'rwc (kg/m^3)'         ,'iwc (kg/m^3)'        ,'rain_gate (kg/m^2/s)',$
    'rain_local (kg/m^2/s)','spectral width (m/s)','u_7r (m/s)'          ,$
    'v_7r (m/s)'           ,'wspd_7r (m/s)'       ,'div_1r (1/s)'        ,$
    'div_3r (1/s)'         ,'div_5r (1/s)'        ,'std_div_1r (1/s)'    ,$
    'std_div_3r (1/s)'     ,'std_div_5r (1/s)'    ,'chi_fit_1r'          ,$
    'chi_fit_3r'           ,'chi_fit_5r'          ,'az_obs_1r (deg)'     ,$
    'az_obs_3r (deg)'      ,'az_obs_5r (deg)'     ,'gapmax_1r (deg)'     ,$
    'gapmax_3r (deg)'      ,'gapmax_5r (deg)'     ,'std_u_7r (m/s)'      ,$
    'std_v_7r (m/s)'       ,'chi_fit_7r']
