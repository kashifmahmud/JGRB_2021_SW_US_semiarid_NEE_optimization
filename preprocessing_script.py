#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 10:00:00 2020

@author: kashifmahmud

Semi-arid site analysis - Script to pre-process the netcdf climate forcing file 
to merge flux data for all US Semi-arid flux sites 
Data - Both forcing netcdf and flux text files are from Natasha

"""

# import necessary libraries
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import datetime
import xarray as xr
import os
import cftime

# ---------------------
# necessary snippets

# snippet to copy a netcdf file to another file
def create_file_from_source(src_file, trg_file):
    src = Dataset(src_file)
    trg = Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]

    # Save the file
    trg.close()


#%%
# ---------------------
# -
# - Set-up
# -
rundir = '/Users/kashifmahmud/IU/Work/All_analyses/semi-arid_sites/' # - CHANGE THIS
site_file = rundir+'siteinfo.txt' # - SET-UP THIS FILE FOR ALL SITES (SITE NAME, START YEAR, END YEAR)
#site_file = rundir+'siteinfo_ses_seg.txt' # - SET-UP THIS FILE FOR ALL SITES (SITE NAME, START YEAR, END YEAR)

# - make output directory(ies)
outdir = rundir+'output/'
if not os.path.isdir(outdir): os.mkdir(outdir)
if not os.path.isdir(outdir+'data/'): os.mkdir(outdir+'data/')

site_names = np.loadtxt(site_file, dtype=str)[:,0]
forcing_starts = np.asarray(np.loadtxt(site_file, dtype=str)[:,1], dtype=int)
forcing_ends = np.asarray(np.loadtxt(site_file, dtype=str)[:,2], dtype=int)
    
    
# - modify the names of input and output flux and forcing files
input_flux_root = rundir + 'data/%s_fluxes.txt'
leap_forcing_nc4_root = rundir + 'data/%s_forcing_leap.nc'
leap_forcing_nc3_root = rundir + 'data/%s_forcing_leap_nc3.nc'
noleap_forcing_root = rundir + 'data/%s_forcing_noleap.nc'

input_forcing_root = rundir + 'output/data/%s_forcing_noleap_nolev.nc'
out_flux_root = rundir + 'output/data/%s_fluxes.csv'
out_forcing_root = rundir + 'output/data/%s_%(start)s-%(end)s_final.nc'
out_cropped_root = rundir + 'output/data/%s_%(start)s-%(end)s_cropped.nc'


#%%
# ---------------------
# -
# - Run the loop for all sites listed in 'siteinfo.txt'
# -
for ns, ss in enumerate(site_names):
    
    #%%
# ---------------------

#    ns = 1
#    ss = site_names[ns]
    # - get Ameriflux data
    input_flux = input_flux_root.replace("%s", ss)
    leap_forcing_nc4 = leap_forcing_nc4_root.replace("%s", ss)
    leap_forcing_nc3 = leap_forcing_nc3_root.replace("%s", ss)
    noleap_forcing = noleap_forcing_root.replace("%s", ss)
    input_forcing = input_forcing_root.replace("%s", ss)
    out_flux = out_flux_root.replace("%s", ss)
    out_forcing = out_forcing_root.replace("%s", ss)
    out_forcing = out_forcing % dict(start = forcing_starts[ns], end = forcing_ends[ns])
    out_cropped = out_cropped_root.replace("%s", ss)

    # - run the terminal command from python to convert netcdf4 to netcdf (classic/nc3)
    os.system('nccopy -k 1 ' + leap_forcing_nc4 + ' ' + leap_forcing_nc3)
    
    # - run the terminal command from python to remove all leap days (Feb. 29) from netcdf3 files
    os.system('cdo delete,month=2,day=29 ' + leap_forcing_nc3 + ' ' + noleap_forcing)
    
    # - run the terminal command from python to remove the lev (heigt) dimension from 
    # the netcdf forcing file using NCO command ncwa
    os.system('ncwa -a lev ' + noleap_forcing + ' -O ' + input_forcing)
    
#    data_flux_v0 = Dataset(input_forcing, 'a')
#    tstep = data_flux_v0.variables['tstep']
#    tstep.calendar = "gregorian"
#    data_flux_v0.close()

    # 140160 / 365/8/48
    
    # ---------------------
    # - load flux text data for US-Whs site
    flux_raw = pd.read_csv(input_flux, sep=" ")
    
    # - subset the NEE, GPP, Reco
    flux = flux_raw[['#Date','NEE','GPP','Reco']]

    # - rename column
    flux = flux.rename(columns = {"#Date":"date"})
    
    # - remove all 29 Febs from the data frame
    flux['date'] = pd.to_datetime(flux['date'])
    flux_noleap = flux[~((flux.date.dt.month == 2) & (flux.date.dt.day == 29))]
    
    # - save the final flux data in csv
    flux_noleap.to_csv(out_flux, index=False)
    
    # - convert the GPP values from float64 to float32 to match with the netcdf file
    flux_noleap = pd.read_csv(out_flux, dtype={col: np.float32 for col in ['NEE','GPP','Reco']})
    flux_noleap['date'] =  pd.to_datetime(flux_noleap['date'], format='%Y-%m-%d')
    
    
    # ---------------------
    #read the forcing netcdf file using xarray to find the starting date, end date and data freequency
    forcing_data = xr.open_dataset(input_forcing)
    if forcing_data['tstep'].dtype == '<M8[ns]': 
        forcing_time = forcing_data.indexes['tstep'] # for US-Fuf
    else: 
        forcing_time = forcing_data.indexes['tstep'].to_datetimeindex() # for other sites
    
    forcing_freq = (forcing_time[1] - forcing_time[0]).total_seconds()
    forcing_freq = str(int(forcing_freq/60))+'T'
    
    # - generate a data frame with same time axis as the flux netcdf data
    timestamp = pd.date_range(start=forcing_time[0], end=forcing_time[-1], freq=forcing_freq)
    date = pd.to_datetime([d.date() for d in timestamp])
    flux_new = pd.DataFrame({'timestamp':timestamp,'date':date})
    
    # - remove all 29 Febs from the data frame
    flux_new = flux_new[~((flux_new.timestamp.dt.month == 2) & (flux_new.timestamp.dt.day == 29))]
    
    # - merge the data frames to format GPP having forcing time steps
#    flux_merge = pd.merge(left=flux_new,right=flux_noleap, how='outer', left_on='date', right_on='date')       
    flux_merge = pd.merge(left=flux_new,right=flux_noleap, how='left', left_on='date', right_on='date')       
    flux_final = flux_merge.replace(-9999.0, np.nan)
    
    
    # ---------------------
    # - make a copy of the netcdf file using snippet (this option is faster)
    create_file_from_source(input_forcing, out_forcing)
    data_flux_v1 = Dataset(out_forcing, 'a')
    
#    tstep = data_flux_v1.variables['tstep']
#    print(tstep[0])
#    print(tstep[tstep.size-1])
#    print(tstep.size)
    
    # - create a new variables GPP for the netcdf flux data
    GPP = data_flux_v1.createVariable("GPP","f4",("tstep","lat","lon"))
    GPP.long_name = "Gross Primary Production"
    GPP.units = "gC/m2/tstep"
    
    # - create a new variables NEE for the netcdf flux data
    NEE = data_flux_v1.createVariable("NEE","f4",("tstep","lat","lon"))
    NEE.long_name = "Net Ecosystem Exchange"
    NEE.units = "gC/m2/tstep"
    
    # - create a new variables Reco for the netcdf flux data
    Reco = data_flux_v1.createVariable("Reco","f4",("tstep","lat","lon"))
    Reco.long_name = "Ecosystem Respiration"
    Reco.units = "gC/m2/tstep"
    
    # - assign the GPP values to the variable in the netcdf file using a loop
    for nt in range(0,len(GPP)):
        GPP[nt] = flux_final['GPP'][nt]
        NEE[nt] = flux_final['NEE'][nt] 
        Reco[nt] = flux_final['Reco'][nt] 
        
#    flux_final['GPP'][0]
#    flux_final['GPP'][140159]
    
    data_flux_v1.close()

#    # - read the netcdf file and change the time calender to 'gregorian' to match 
#    # with datetime format for cropping in next step
#    data_flux_v2 = Dataset(out_forcing, 'a')
#    tstep = data_flux_v2.variables['tstep']
#    tstep.calendar = "gregorian"
#    data_flux_v2.close()
#
#    gpp = data_flux_v2.variables['GPP']
#    print(gpp[0])
#    print(gpp[gpp.size-1])
#    print(gpp.size)
    

    ## ---------------------
    # - crop the data for the available flux data using xarray
    # - first crop the flux data frame for missing values
    flux_crop = flux_final.dropna(how='any') 
    
    #    #GPP and NEE time series plot to check the data
    #    import matplotlib.pyplot as plt
    #    plt.plot(flux_crop.date,flux_crop.GPP)
    #    plt.plot(flux_crop.date,flux_crop.NEE)
    #    plt.plot(flux_crop.date,flux_crop.Reco)

    #find the starting time point of the flux data
    if flux_crop.timestamp.dt.month.iloc[0] != 1 or flux_crop.timestamp.dt.day.iloc[0] != 1:
        cropping_starts = datetime.datetime((flux_crop.timestamp.dt.year.iloc[0] + 1), 1, 1)
#        cropping_starts = cropping_starts.strftime("%Y-%m-%d %H:%M:%S")
        
    else:
#        cropping_starts = flux_crop.timestamp.iloc[0]
        cropping_starts = flux_crop.timestamp.iloc[0].to_pydatetime()
    
    #find the final time point of the flux data
    if flux_crop.timestamp.dt.month.iloc[-1] != 12 or flux_crop.timestamp.dt.day.iloc[-1] != 31:
        cropping_ends = datetime.datetime((flux_crop.timestamp.dt.year.iloc[-1] - 1), 12, 31, 23, 30)
#        cropping_ends = cropping_ends.strftime("%Y-%m-%d %H:%M:%S")
    else:
#        cropping_ends = flux_crop.timestamp.iloc[-1]
        cropping_ends = flux_crop.timestamp.iloc[-1].to_pydatetime()
    
    # - read the final netcdf file
    data = xr.open_dataset(out_forcing)
#    data['tstep']
    tstep = pd.DataFrame(data['tstep'].values[:])
    
    if data['tstep'].dtype == 'object':
    # convert both start and end time stamps to noleap datetime to match with the netcdf time stamp
        cropping_starts = cftime.DatetimeNoLeap(cropping_starts.year,cropping_starts.month,cropping_starts.day,0,0,0)
        cropping_ends = cftime.DatetimeNoLeap(cropping_ends.year,cropping_ends.month,cropping_ends.day,cropping_ends.hour,cropping_ends.minute,cropping_ends.second) 
    
    #crop the netcdf file
    data_crop = data.sel(tstep=slice(cropping_starts, cropping_ends))
#    data_crop['tstep']
#    tstep = pd.DataFrame(data_crop['tstep'].values[:])
    
    
    #save the cropped data as netcdf file
    out_cropped = out_cropped % dict(start = cropping_starts.strftime("%Y"), end = cropping_ends.strftime("%Y"))
    data_crop.to_netcdf(out_cropped)
    
    # - read the netcdf file and change back the time calender to 'noleap'
    data_flux_cropped = Dataset(out_cropped, 'a')
    tstep = data_flux_cropped.variables['tstep']
    tstep[:] = tstep[:] - tstep[0]
#    tstep[:] = tstep[:] - 3600*24*365*(cropping_starts.year - flux_final.date[0].year)
#    print(tstep[tstep.size-1])
#    print(tstep.size)
    tstep.calendar = "noleap"
    tstep.units = 'seconds since ' + cropping_starts.strftime("%Y-%m-%d")
    data_flux_cropped.close()

#    ss
#
    data_flux_cropped = Dataset(out_cropped, 'a')
    tstep = data_flux_cropped.variables['tstep']
    del tstep._FillValue
    gpp = data_flux_cropped.variables['GPP']
    del gpp._FillValue
    nee = data_flux_cropped.variables['NEE']
    del nee._FillValue
    reco = data_flux_cropped.variables['Reco']
    del reco._FillValue
    data_flux_cropped.close()
    
#    np.argwhere(np.isnan(gpp[:]))
#    gpp[5215]
#    gpp[5216]

#data_flux_v0 = Dataset(input_forcing, 'a')
#tstep = data_flux_v0.variables['tstep']
#data_flux_v0.close()
#
#252459000 - 3600*24*365*8
#252459000 - (365*6 + 366*2)*24*3600
#
#3.162240e+07 - 3600*24*365
#86400/3600
#
#(105120)/365/6
#365*6*48
