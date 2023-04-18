#################################################################################################
#		
# Title:	Anthropogenic climate change impacts exacerbate summer forest fires in California
# 				
# Authors: 	Marco Turco, University of Murcia (marco.turco@um.es) 
#               John T. Abatzoglou, Sixto Herrera, Yizhou Zhuang, Sonia Jerez, Donald D. Lucas, 
#               Amir AghaKouchak6 and Ivana Cvijanovic
#
#################################################################################################

#################################################################################################
# A. General instructions 
#################################################################################################

This project is designed to be executed with shell scripts, R and Matlab codes. 
Execute script files in the order they are listed.

Data sources:

- nClimGrid tmax and prec
https://www.ncei.noaa.gov/data/nclimgrid-monthly/access/

- PRISM climate data tmax, tmin, tdmean 
https://www.prism.oregonstate.edu/

- FRAP fires data
https://frap.fire.ca.gov/media/ly2jyr4j/fire21_2.zip

- MTBS fires data
https://mtbs.gov/direct-download


Notes regarding reproducibility:

Script files starting with "1_" in their name are for data preprocessing. 
Most of these script files will NOT run because we do not include the 
raw data because files are simply too large to conveniently share. 
We suggest you run script starting with "2_" which directly reproduces the 
results in the paper. 

If you have any questions or wish to express any comment to the authors, please 
contact Dr. Marco Turco at the email indicated above.


#################################################################################################
# B. Description of script files
#################################################################################################

Scripts for data preparation and cleaning.

- 1_1_prepare_nclimgrid.sh
Get nclimgrid dataset and select time-space domains.

- 1_2_prepare_nclimgrid.R
Prepare region-level nclimgrid dataset.

- 1_3_prepare_prism.R
Get PRISM data and calculate VPD.

- 1_4_prepare_fire_data.R
Prepare mtbs fires data and frap data without the forest filter.

- 1_5_prepare_frap_data.m
Prepare frap data and compare it with other datasets. Plot figure S1.

- 1_6_prepare_gcm_attribution.R
Prepare regional averages of GCM data.

- 1_7_prepare_gcm_attribution.m
Prepare GCM data for attribution analyses. Plot figure S3.

- 1_8_prepare_gcm_scenarios.R
Prepare GCM data for scenarios analyses.

- 1_9_calculate_future_scenarios.m
Computes future impacts of anthropogenic climate change.

- 1_10_calculate_future_scenarios_tp.m
Computes future impacts of anthropogenic climate change with the model that includes the precipitation effect.

Scripts to reproduce the results in the paper.

- 2_1_climate_fire_model.R
Build and test the best climate-fire model through out-of-sample calibration.

- 2_2_climate_fire_model.m
Validate the best climate-fire model through out-of-sample tenfold validation. Plot figure 1 and S2.

- 2_3_attribution.m
Compute impacts of anthropogenic climate change. Plot figure 2 and S4.

- 2_4_attribution_tp.m
Compute anthropogenic climate change impacts with the model that includes the precipitation effect. 
Plot figure S5.

- 2_5_plot_future_scenarios.m
Compute future BA simulations. To run this script you have to download the data from: 
https://univmurcia-my.sharepoint.com/:f:/r/personal/marco_turco_um_es/Documents/gcms_pnas_2022?csf=1&web=1&e=ak4Qql

- 2_5_plot_future_scenarios_tp.m
Compute future BA simulations with the model that includes the precipitation effect. Plot figure 3.
To run this script you have to download the data from: 
https://univmurcia-my.sharepoint.com/:f:/r/personal/marco_turco_um_es/Documents/gcms_pnas_2022?csf=1&web=1&e=ak4Qql


In the folder "misc" there are additional script files useful for the execution of the main script files described above.