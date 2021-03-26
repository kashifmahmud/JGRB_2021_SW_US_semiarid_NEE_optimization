# ORCHIDEE parameter optimization for SW US semiarid sites
Steps to follow:
1. The meteorological forcing data and eddy covariance measurements of net surface energy and carbon exchanges at 30-minutes intervals needed to run the model optimization are freely available from figshare (https://doi.org/10.6084/m9.figshare.14327489.v3). Download all the data files and save in the "data" folder. 
2. Run the python script "all_semiarid.py" to pre-process the netcdf climate forcing file to merge flux data for all US semiarid flux sites.
3. Run the ORCHIDAS simulation with the processed netcdf files for all sites. The model output files can be downloaded from the figshare.
4. Save the model output files in the "output" folder. 
5. Finally run the python script "plot_outputs_semiarid_final.py" to plot the figures and tables presented in the manuscript using ORCHIDAS optimization outputs.
