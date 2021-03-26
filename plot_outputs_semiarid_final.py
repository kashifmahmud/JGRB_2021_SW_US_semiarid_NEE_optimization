#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 14 14:20:51 2020

@author: kashifmahmud

Plot the output files from Orchidas optimizations for all semi-arid sites

"""

# import necessary libraries
import pandas as pd
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import linregress
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import statsmodels.formula.api as smf
from matplotlib.patches import Rectangle

#from math import pi
#import matplotlib.lines as mlines            
#from netCDF4 import Dataset
#import xarray as xr

## ---------------------
# necessary snippets
# snippet to calculate Mean square deviation (MSD) based on Kobayoshi and Salam (2000)
def decomp_mse_kobayashi(x, y):
        """
        Stats Kobayashi
        Author: Cedric Bacour
        """
        X = x.ravel(); Y = y.ravel()
        ix=np.where(~np.isnan(X))[0];iy=np.where(~np.isnan(Y))[0]
        ii = list(set(ix).intersection(iy))
        X=X[ii];Y=Y[ii]
        sX = np.nanstd(X)
        sY = np.nanstd(Y)
        r = corrcoef_func(X,Y)
        bias = np.nanmean(X-Y)**2
        variance = (sX-sY)**2
        phase = 2*(sX*sY)*(1-r)
        mse = (calc_rmse(X,Y))**2
        return mse,bias,variance,phase,r
    
# snippet to calculate correlation coefficient
def corrcoef_func(x, y, npixMin = None):
        """
		Correlation function
		Author: Cedric Bacour
		"""
        iOK = np.where((~np.isnan(x.ravel())==True) & (~np.isnan(y.ravel())==True))
        if npixMin is not None:
            if len(iOK[0]) < npixMin: return np.nan
        return np.corrcoef(x.ravel()[iOK],y.ravel()[iOK])[1,0]
	
# snippet to calculate RMSE
def calc_rmse(x, y):
        """
		RMSE function
		Author: Cedric Bacour
		"""
        iOK = np.where((~np.isnan(x)==True) & (~np.isnan(y)==True))
        n = len(iOK[0])
        rmse = np.sqrt( np.sum((x[iOK]-y[iOK])**2)/n )
        return rmse

##snippet to draw spider plot
#def make_spider( row, title, color, optim, param_title, i, fig_title, nf):
#    # number of variable
#    categories=list(monthly_ST)[1:]
#    N = len(categories)
#    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
#    angles = [n / float(N) * 2 * pi for n in range(N)]
#    angles += angles[:1]
#    # Initialise the spider plot
#    ax = plt.subplot(3,4,row+1, polar=True, )
#    # If you want the first axis to be on top:
#    ax.set_theta_offset(pi / 2)
#    ax.set_theta_direction(-1)
#    # Draw one axe per variable + add labels labels yet
#    plt.xticks(angles[:-1], categories, color='grey', size=8)
#    # Draw ylabels
#    ax.set_rlabel_position(0)
#    plt.yticks([0.5,1], ["0.5","1.0"], color="grey", size=8)
#    plt.ylim(0,1)
#    # Ind1
#    values=monthly_ST.loc[row].drop('Site').values.flatten().tolist()
#    values += values[:1]
#    ax.plot(angles, values, color=color, linewidth=1, linestyle='solid')
#    ax.fill(angles, values, color=color, alpha=0)
#    #Add a title
#    plt.title(title, size=11, color='black', y=1.1)
#    plt.suptitle(fig_title[nf]+' Taylor skill score (ST) ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i]) 
       
             
                
#%%
## ---------------------
## ---------------------
# Set-up
rundir = '/Users/kashifmahmud/IU/Work/All_analyses/SW-US_semiarid/output/for_plots/' # - CHANGE THIS
site_file = rundir+'siteinfo.txt' # - SET-UP THIS FILE FOR ALL SITES (SITE NAME, START YEAR, END YEAR)
param_file = rundir+'param_list.xlsx'
flux_names = ['NEE','GPP','Reco']

# - make output directory(ies)
outdir = rundir+'figures_tables/'
if not os.path.isdir(outdir): os.mkdir(outdir)
if not os.path.isdir(outdir + 'Flux_timeseries/'): os.mkdir(outdir + 'Flux_timeseries/')
if not os.path.isdir(outdir + 'Flux_monthly/'): os.mkdir(outdir + 'Flux_monthly/')
if not os.path.isdir(outdir + 'Flux_annual_anomaly/'): os.mkdir(outdir + 'Flux_annual_anomaly/')
if not os.path.isdir(outdir + 'RMSE/'): os.mkdir(outdir + 'RMSE/')
if not os.path.isdir(outdir + 'Mean_monthly_fluxes/'): os.mkdir(outdir + 'Mean_monthly_fluxes/')
if not os.path.isdir(outdir + 'MSD/'): os.mkdir(outdir + 'MSD/')
if not os.path.isdir(outdir + 'Parameters/'): os.mkdir(outdir + 'Parameters/')
#if not os.path.isdir(outdir + 'Taylor_skill_score/'): os.mkdir(outdir + 'Taylor_skill_score/')
if not os.path.isdir(outdir + 'Flux_corr/'): os.mkdir(outdir + 'Flux_corr/')
if not os.path.isdir(outdir + 'Flux_annual_anomaly/box_whisker_plots/'): os.mkdir(outdir + 'Flux_annual_anomaly/box_whisker_plots/')
if not os.path.isdir(outdir + 'Scatter_plots/'): os.mkdir(outdir + 'Scatter_plots/')
if not os.path.isdir(outdir + 'Manuscript_figures/'): os.mkdir(outdir + 'Manuscript_figures/')
 
# define the optimization simulations for plotting
#optim_name = ['gpp_reco_iter3840_paramA','nee_iter3840_paramA']
#param_title = ['All parameters','All parameters']
#optim = ['(GPP + $R_{eco}$)','NEE']
#optim_name = ['gpp_reco_iter3840_paramA','gpp_reco_iter2400_paramB','gpp_reco_iter3360_paramC','gpp_reco_iter1680_paramD',
#              'gpp_reco_iter2160_paramE','gpp_reco_iter720_paramF','nee_iter3840_paramA','nee_iter2400_paramB',
#              'nee_iter3360_paramC','nee_iter1680_paramD','nee_iter2160_paramE','nee_iter720_paramF',
#              'gpp_iter3360_paramC','gpp_iter1680_paramD','gpp_iter2160_paramE',
#              'reco_iter3840_paramA','reco_iter3360_paramC']
#param_title = ['All parameters','Photosynthesis and post-GPP parameters','Photosynthesis and phenology parameters',
#               'Photosynthesis parameters','Phenology parameters','Post-GPP parameters',
#               'All parameters','Photosynthesis and post-GPP parameters','Photosynthesis and phenology parameters',
#               'Photosynthesis parameters','Phenology parameters','Post-GPP parameters',
#               'Photosynthesis and phenology parameters','Photosynthesis parameters','Phenology parameters',
#               'All parameters','Photosynthesis and phenology parameters']
#optim = ['(GPP + $R_{eco}$)','(GPP + $R_{eco}$)','(GPP + $R_{eco}$)','(GPP + $R_{eco}$)','(GPP + $R_{eco}$)','(GPP + $R_{eco}$)',
#         'NEE','NEE','NEE','NEE','NEE','NEE','GPP','GPP','GPP','$R_{eco}$','$R_{eco}$']
#param_group = {'param':['Prior','P1','P2','P3','P4','P5','P6','P1','P2','P3','P4','P5','P6','P3','P4','P5','P1','P3'], 'optim':[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]}

optim_name = ['nee_iter3840_paramP1','nee_iter3360_paramP2','nee_iter2880_paramP3',
              'nee_iter2400_paramP4','nee_iter2160_paramP5','nee_iter1680_paramP6','nee_iter720_paramP7']
param_title = ['All parameters','Phenology and Photosynthesis parameters','Phenology and Post-GPP parameters',
               'Photosynthesis and Post-GPP parameters','Phenology parameters','Photosynthesis parameters','Post-GPP parameters']
optim = ['NEE','NEE','NEE','NEE','NEE','NEE','NEE']
param_group = {'param':['Prior','P1','P2','P3','P4','P5','P6','P7'], 
               'optim':[0,1,2,3,4,5,6,7],
               'opt':['Prior','Opt1','Opt1','Opt1','Opt1','Opt1','Opt1','Opt1']}

fig_title = ['NEE','GPP','$R_{eco}$']
subplot_title = ['(a)','(b)','(c)','(d)','(e)','(f)']
param_plot = True
rmse_plot = True
corr_plot = True
ms_plot = True
scatter_plot = True

# - set figure styles
sns.set_style("ticks")
plt.style.use('seaborn-white')
pal = dict(Obs="gray", Prior="green", Post="red")
col_wrap = 3
    
df_rmse_all = df_param_all = flux_corr_all = flux_slope_all = flux_msd_all = flux_msd_all_month = flux_msd_all_annual = pd.DataFrame([])
df_annual_all = df_month_mean_all = pd.DataFrame([])

#%%
for i in range(len(optim_name)):
#    i=0
    # - modify the names of input and output flux and forcing files
    flux_root = rundir + '%s/%s_' + optim_name[i] + '/output.nc'
    #flux_root = rundir + '%s/%s_nee_iter2250_paramA/output.nc'
    site_names = np.loadtxt(site_file, dtype=str)[:,0]
    index = np.argwhere(site_names=='NA')
    site_names = np.delete(site_names, index)
    df_site_names = pd.DataFrame(site_names)
    df_site_names.columns = ['Site']
                       
    veg_type = np.loadtxt(site_file, dtype=str)[:,1]
    veg_type = np.delete(veg_type, index)
    
    var_name_root = 'data_site0_var%s'
    
    #creating the fluxes and RMSE data frames
    ## ---------------------
    ## ---------------------
    # initialize the data frames and arrays
    df_flux = df_rmse = df_rmse_site = df_anomaly = anomaly_site_all = df_month_mean = df_error = pd.DataFrame([])
    df_param = df_param_sub = df_msd = df_msd_month = df_msd_annual = df_month = df_annual = pd.DataFrame([])
    plot_nee_text = plot_gpp_text = plot_reco_text = []
    
    for ns, ss in enumerate(site_names):  
    #    ns = 1
        ss = site_names[ns]
        vt = veg_type[ns]
        flux = flux_root.replace("%s", ss)
        
        ###### GPP and NEE time series plot
        # - load output nc file
        data = netCDF4.Dataset(flux, 'r')
        
        ## ---------------------
        # get the time steps for each site
        NEE = data.variables['data_site0_var0']
        time = netCDF4.num2date(float(NEE.getncattr("time_offset")) + np.arange(len(NEE[0,:])) * float(NEE.getncattr("time_step")), units = NEE.getncattr("time_units"))
        time = pd.to_datetime(time.astype('str'))
        years_df = pd.DatetimeIndex(time)
#        years_df = time.to_datetimeindex()
        years = years_df.year
        years = years.unique()
        
        
        ## ---------------------
        # get the rmse for each site
        rmse = data.variables['RMSE'][[0,-1],:]
        rmse_site = np.concatenate((np.array(rmse[0,:][:]), np.array(rmse[1,:][:])))
        nee_text = ['%s' % ss + ' (' + '%s' % vt + ') - ' + 'RMSE: Prior = %0.3f' % rmse_site[0][0] + ', ' + 'Post = %0.3f' % rmse_site[1][0]]
        gpp_text = ['%s' % ss + ' (' + '%s' % vt + ') - ' + 'RMSE: Prior = %0.3f' % rmse_site[0][1] + ', ' + 'Post = %0.3f' % rmse_site[1][1]]
        reco_text = ['%s' % ss + ' (' + '%s' % vt + ') - ' + 'RMSE: Prior = %0.3f' % rmse_site[0][2] + ', ' + 'Post = %0.3f' % rmse_site[1][2]]
        
        # appending the text array with all site RMSE
        plot_nee_text = np.concatenate((plot_nee_text,nee_text))
        plot_gpp_text = np.concatenate((plot_gpp_text,gpp_text))
        plot_reco_text = np.concatenate((plot_reco_text,reco_text))    
        
        #    rmse_site = pd.DataFrame(rmse[0,:][:])
        df_rmse_site = pd.DataFrame(rmse[0,:][:])
        df_rmse_site = df_rmse_site.append(pd.DataFrame(rmse[1,:][:]))
        df_rmse_site.loc[2] = df_rmse_site.iloc[1,:] / df_rmse_site.iloc[0,:]
        df_rmse_site.columns = ['NEE','GPP','Reco']
        df_rmse_site['Site'] = ss
        df_rmse_site['Datatype'] = ['Prior','Post','Ratio']
        df_rmse_melt = pd.melt(df_rmse_site, id_vars=('Site','Datatype'), var_name="Flux", value_name="RMSE")
            
        # appending the data frame with all site RMSE
        df_rmse = df_rmse.append(df_rmse_melt)
        
        
        ## ---------------------
        # get the parameters for each site
        param = pd.DataFrame(data.variables['param_id'][:])
        param.columns = ['Param_ID']
        param['Param_post'] = data.variables['param'][[-1],:].transpose()
        param['Param_default'] = data.variables['param_default'][:]
        param['Param_min'] = data.variables['param_min'][:]
        param['Param_max'] = data.variables['param_max'][:]
            
        # take the diagonal values of 'param_bpost' for posterior parameter uncertainty
        param['Param_error'] = (np.diag(pd.DataFrame(data.variables['param_bpost'][:]))) ** 0.5
        param['Param_error_max'] = param['Param_post'] + param['Param_error']
        param['Param_error_min'] = param['Param_post'] - param['Param_error']
        
#        #subset few parameters to display in figures 
#        if 'paramA' in optim_name[i]: 
#            param_sub = param[param['Param_ID'].str.contains(
#                    'VCMAX25|HYDROL_HUMCSTE|TPHOTO_MAX|TPHOTO_MIN|SLA|LAI_MAX|LEAFAGECRIT|TAU_LEAFINIT|LAI_MAX_TO_HAPPY|MIN_GROWTHINIT_TIME|KSOILC|SOIL_Q10|FRAC_GROWTHRESP')]
        
        param['Param_prior_error'] = data.variables['param_error'][:] # parameter prior uncertainty values
        param['deviation'] = (param['Param_post'] - param['Param_default']) / (param['Param_max'] - param['Param_min']) 
        param['uncer_reduc'] = 1 - (param['Param_error'] / param['Param_prior_error'])
        
        param.index = np.arange(len(param))
        param['Param_name'] = param['Param_ID'].str.split('__').str[0]
        param_melt = pd.melt(param, id_vars=('Param_ID','Param_name'), var_name="Paramtype", value_name="Paramvalue")
        param_melt['Site'] = ss
        # appending the data frame with all site parameters
        df_param = df_param.append(param_melt)
        
        
#        param_sub.index = np.arange(len(param_sub))
#        param_sub['Param_name'] = param_sub['Param_ID'].str.split('__').str[0]
#        param_sub_melt = pd.melt(param_sub, id_vars=('Param_ID','Param_name'), var_name="Paramtype", value_name="Paramvalue")
#        param_sub_melt['Site'] = ss
#        # appending the data frame with all site parameters
#        df_param_sub = df_param_sub.append(param_sub_melt)


        ## ---------------------
        # extract all 3 flux data into dataframe for each site
        for nf, sf in enumerate(flux_names):     
            # nf = 0
            sf = flux_names[nf]
            
            # Read the flux variables
            var_name = var_name_root.replace("%s", str(nf))
            dim = data.variables[var_name][:].shape
            if dim[0] == 4:
                flux_site = data.variables[var_name][[0,1,3], :]
            elif dim[0] == 3: # US-Fuf has one less dimension in the data
                flux_site = data.variables[var_name][[0,1,2], :]
            
            ## ---------------------
            # Read the errors
            error_site = pd.DataFrame(data.variables['site_error'][:])
            error_site.columns = ['NEE','GPP','Reco']
            error_site['Site'] = ss
            # append the dataframe with new flux data
            df_error = df_error.append(error_site)
            
        
            ## ---------------------
            #create annual flux anomaly data frame 
            data_annual = np.sum(np.reshape(flux_site[0,:], (len(years),365)), axis = 1)
            data_anomaly = data_annual - np.mean(data_annual)
                
            prior_annual = np.sum(np.reshape(flux_site[1,:], (len(years),365)), axis = 1)
            prior_anomaly = prior_annual - np.mean(prior_annual)
                
            post_annual = np.sum(np.reshape(flux_site[2,:], (len(years),365)), axis = 1)
            post_anomaly = post_annual - np.mean(post_annual)
    
            anomaly_site = pd.DataFrame({'Obs': data_anomaly, 'Prior': prior_anomaly, 'Post': post_anomaly})
            anomaly_site['Year'] = years
            anomaly_site['Site'] = ss
            anomaly_site['Flux'] = sf
            anomaly_site_all = anomaly_site_all.append(anomaly_site)
            
            anomaly_melt = pd.melt(anomaly_site, id_vars=('Year','Site','Flux'), var_name="Datatype", value_name="Anomaly")
            # append the dataframe with new flux data
            df_anomaly = df_anomaly.append(anomaly_melt)
            
            
            ## ---------------------
            #create daily flux data frame 
            flux_site = pd.DataFrame(flux_site)
            flux_site = flux_site.transpose()
            flux_site.columns = ['Obs','Prior','Post']
            flux_site['Date'] = time
            flux_site['Date'] = flux_site['Date'].dt.date 
            flux_melt = pd.melt(flux_site, id_vars='Date')
            flux_melt['Site'] = ss
            flux_melt['Flux'] = sf
            
            # append the dataframe with new flux data
            df_flux = df_flux.append(flux_melt)
    
            
            ## ---------------------
            #create monthly flux data frame 
            flux_site['Date'] = pd.to_datetime(flux_site['Date']) # convert to datetime
            data_month = flux_site.set_index('Date').groupby(pd.Grouper(freq='M'))['Obs','Prior','Post'].sum().reset_index()
            data_month['Month'] = data_month['Date'].dt.month
            data_month_mean = data_month.groupby('Month', as_index=False).mean()
            data_month_mean = pd.melt(data_month_mean, id_vars='Month', var_name="Datatype", value_name="Monthly_mean")
            data_month_mean['Site'] = ss
            data_month_mean['Flux'] = sf
            # append the dataframe with new flux data
            df_month_mean = df_month_mean.append(data_month_mean)
    
            #create total monthly flux data frame
            data_month['Site'] = ss
            data_month['Flux'] = sf
            # append the dataframe with new flux data
            df_month = df_month.append(data_month)
#            df_month_mean_unstack = df_month.set_index(['Site','Flux','Month','Datatype'])['Monthly_mean'].unstack().reset_index()

            ## ---------------------
            #create annual flux data frame 
            data_annual = flux_site.set_index('Date').groupby(pd.Grouper(freq='Y'))['Obs','Prior','Post'].sum().reset_index()
            data_annual['Year'] = data_annual['Date'].dt.year
            data_annual['Site'] = ss
            data_annual['Flux'] = sf
            # append the dataframe with new flux data
            df_annual = df_annual.append(data_annual)
    
            
            ## ---------------------
            # Calculate Mean square deviation (MSD) for daily data based on Kobayoshi and Salam (2000)
            msd_site = pd.DataFrame(decomp_mse_kobayashi(flux_site.Obs, flux_site.Post))
            msd_site[1] = pd.DataFrame(decomp_mse_kobayashi(flux_site.Obs, flux_site.Prior))
            msd_site = msd_site.transpose()
            msd_site.columns = ['mse','bias','variance','phase','r']            
            msd_site['Sim'] = ['post','prior']
            msd_site['Site'] = ss
            msd_site['Flux'] = sf
            # append the dataframe with new statistics data
            df_msd = df_msd.append(msd_site)
            
            ## ---------------------
            # Calculate Mean square deviation (MSD) for monthly data based on Kobayoshi and Salam (2000)
            msd_site_month = pd.DataFrame(decomp_mse_kobayashi(data_month.Obs, data_month.Post))
            msd_site_month[1] = pd.DataFrame(decomp_mse_kobayashi(data_month.Obs, data_month.Prior))
            msd_site_month = msd_site_month.transpose()
            msd_site_month.columns = ['mse','bias','variance','phase','r']            
            msd_site_month['Sim'] = ['post','prior']
            msd_site_month['Site'] = ss
            msd_site_month['Flux'] = sf
            # append the dataframe with new statistics data
            df_msd_month = df_msd_month.append(msd_site_month)
            
            ## ---------------------
            # Calculate Mean square deviation (MSD) for annual data based on Kobayoshi and Salam (2000)
            msd_site_annual = pd.DataFrame(decomp_mse_kobayashi(data_annual.Obs, data_annual.Post))
            msd_site_annual[1] = pd.DataFrame(decomp_mse_kobayashi(data_annual.Obs, data_annual.Prior))
            msd_site_annual = msd_site_annual.transpose()
            msd_site_annual.columns = ['mse','bias','variance','phase','r']            
            msd_site_annual['Sim'] = ['post','prior']
            msd_site_annual['Site'] = ss
            msd_site_annual['Flux'] = sf
            # append the dataframe with new statistics data
            df_msd_annual = df_msd_annual.append(msd_site_annual)
            
    # appending the RMSE data frame with all optimization and site RMSE
    df_rmse['optim'] = i+1
    df_rmse['optim_name'] = [j.lower().replace(',', '').split('_iter') for j in optim_name][i][0]
    df_rmse_all = df_rmse_all.append(df_rmse)
    
    # appending the parameter deviation data frame with all optimizations and sites
    df_param['optim'] = i+1
    df_param['optim_name'] = [j.lower().replace(',', '').split('_iter') for j in optim_name][i][0]
    df_param_all = df_param_all.append(df_param)
    
    # appending the MSD dataframe for daily data with all optimizations
    df_msd['optim'] = i+1
    df_msd['optim_name'] = [j.lower().replace(',', '').split('_iter') for j in optim_name][i][0]
    flux_msd_all = flux_msd_all.append(df_msd)
    
    # appending the MSD dataframe for monthly data with all optimizations
    df_msd_month['optim'] = i+1
    df_msd_month['optim_name'] = [j.lower().replace(',', '').split('_iter') for j in optim_name][i][0]
    flux_msd_all_month = flux_msd_all_month.append(df_msd_month)
    
    # appending the MSD dataframe for annual data with all optimizations
    df_msd_annual['optim'] = i+1
    df_msd_annual['optim_name'] = [j.lower().replace(',', '').split('_iter') for j in optim_name][i][0]
    flux_msd_all_annual = flux_msd_all_annual.append(df_msd_annual)
    
    # appending the annual flux dataframe with all optimizations
    df_annual_flux = df_annual.copy()
    df_annual_flux['optim'] = i+1
    df_annual_flux['optim_name'] = [j.lower().replace(',', '').split('_iter') for j in optim_name][i][0]
    df_annual_all = df_annual_all.append(df_annual_flux)
    
    # appending the mean monthly flux data frame with all optimizations
    df_month_mean['optim'] = i+1
    df_month_mean['optim_name'] = [j.lower().replace(',', '').split('_iter') for j in optim_name][i][0]
    df_month_mean_all = df_month_mean_all.append(df_month_mean)
    
    ##%% plotting the parameters with range, error bars
    ## ---------------------
    ## ---------------------
    if param_plot == True:
#    if i == 0 or i == 1:
        # choose the parameters for PFT=11 (C4G=85%) for the US-Wjs site
        if 'US-Wjs' in df_param.Site.unique():
            param_wjs = df_param[df_param.Site == 'US-Wjs']
            param_wjs = param_wjs[~param_wjs.Param_ID.str.contains('__04')] # remove PFT=04 (15%) for US-Wjs
            df_param = df_param[~df_param.Site.str.contains('Wjs')]
            df_param = pd.concat([df_param,param_wjs], ignore_index=True)
        
        # choose the parameters for PFT=11 (C4G=44%) for the US-SRG site
        if 'US-SRG' in df_param.Site.unique():
            param_SRG = df_param[df_param.Site == 'US-SRG']
            param_SRG = param_SRG[~param_SRG.Param_ID.str.contains('__06')] # remove PFT=06 (11%) for US-SRG
            df_param = df_param[~df_param.Site.str.contains('SRG')]
            df_param = pd.concat([df_param,param_SRG], ignore_index=True)    
        
        # choose the parameters for PFT=06 (TeBD=35%) for the US-SRM site
        if 'US-SRM' in df_param.Site.unique():
            param_SRM = df_param[df_param.Site == 'US-SRM']
            param_SRM = param_SRM[~param_SRM.Param_ID.str.contains('__11')] # remove PFT=11 (15%) for US-SRM
            df_param = df_param[~df_param.Site.str.contains('SRM')]
            df_param = pd.concat([df_param,param_SRM], ignore_index=True)
        
        # choose the parameters for PFT=05 (TeBE=55%) for the US-Ses site
        if 'US-Ses' in df_param.Site.unique():
            param_Ses = df_param[df_param.Site == 'US-Ses']
            param_Ses = param_Ses[~param_Ses.Param_ID.str.contains('__11')] # remove PFT=11 (25%) for US-Ses
            df_param = df_param[~df_param.Site.str.contains('Ses')]
            df_param = pd.concat([df_param,param_Ses], ignore_index=True)
        
        # choose the parameters for PFT=04 (TeNE=60%) for the US-SRM site
        if 'US-Mpj' in df_param.Site.unique():
            param_Mpj = df_param[df_param.Site == 'US-Mpj']
            param_Mpj = param_Mpj[~param_Mpj.Param_ID.str.contains('__11')] # remove PFT=11 (20%) for US-Mpj
            df_param = df_param[~df_param.Site.str.contains('Mpj')]
            df_param = pd.concat([df_param,param_Mpj], ignore_index=True)
        
        #subset the dataframe
        df_param_plot = df_param[df_param['Paramtype'].isin(['Param_post','Param_default','Param_min','Param_max'])]
        
        # read param list to add the param types
        param_list = pd.read_excel(param_file)
        df_optim = pd.DataFrame({'Site':np.repeat(site_names,len(param_list))})
        param_list = param_list.append([param_list]*(len(site_names)-1),ignore_index=True)
        param_list = pd.concat([param_list,df_optim], axis=1)
        
        # merge the data frames to have param type
        df_param_plot = pd.merge(df_param_plot, param_list, how='right', on=['Param_name','Site'])
        df_param_plot = df_param_plot.sort_values(['Param_type','Param_name'], ascending=True)
    
        #plot all parameters
        sns.set(font_scale=0.8, style='white')
        g = sns.relplot(x='Site', y='Paramvalue', hue='Paramtype', style='Paramtype', 
                size = 'Paramtype', col='Param_name',
                data=df_param_plot, palette=['purple','grey','orange','orange'], 
                col_wrap=11, height=10, aspect=1,
                markers=["o","X",'v','^'], sizes=(10,30),
                facet_kws=dict(sharey=False),
                size_order = ['Param_post','Param_default','Param_min','Param_max'])  

        plt.xlim(-0.5, 12)
        plt.subplots_adjust(hspace=0.3, wspace=0.5)
        [plt.setp(ax.texts, text="") for ax in g.axes.flat]
        g.set_titles(col_template = '{col_name}', size=8.5)
        g.set(xlabel='', ylabel='')
        plt.subplots_adjust(top=0.97,left=0.035,bottom=0.07,right=0.99)
        g._legend.set_bbox_to_anchor([0.99,0.08])
        # replace labels
        new_labels = ['','Posterior parameters','Prior parameters','Parameter bounds','Parameter bounds']
        for t, l in zip(g._legend.texts, new_labels): t.set_text(l)
            
#        g.fig.suptitle('Parameters ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i])   
        g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='Parameter values', size=12, rotation=90)
            
        for nx, ax in enumerate(g.axes.flat):
#            nx = 0
#            ax = g.axes.flat[nx]
            plt.setp(ax.get_xticklabels(), rotation=90, fontsize=7)
                
        g.fig.savefig(outdir + 'Parameters/' + 'Parameters_all_optim_' + optim_name[i] + '_Csink.png', dpi = 1000)
        plt.close()
        
        
#        #plot a subset of parameters
#    if param_plot == True and 'paramA' in optim_name[i]:
#        if 'US-Wjs' in df_param_sub.Site.unique():
#            param_wjs = df_param_sub[df_param_sub.Site == 'US-Wjs']
#            param_wjs = param_wjs[~param_wjs.Param_ID.str.contains('__04')] # remove PFT=04 (15%) for US-Wjs
#            df_param_sub = df_param_sub[~df_param_sub.Site.str.contains('Wjs')]
#            df_param_sub = pd.concat([df_param_sub,param_wjs], ignore_index=True)
#        
#        # choose the parameters for PFT=11 (C4G=44%) for the US-SRG site
#        if 'US-SRG' in df_param_sub.Site.unique():
#            param_SRG = df_param_sub[df_param_sub.Site == 'US-SRG']
#            param_SRG = param_SRG[~param_SRG.Param_ID.str.contains('__06')] # remove PFT=06 (11%) for US-SRG
#            df_param_sub = df_param_sub[~df_param_sub.Site.str.contains('SRG')]
#            df_param_sub = pd.concat([df_param_sub,param_SRG], ignore_index=True)  
#        
#        # choose the parameters for PFT=06 (TeBD=35%) for the US-SRM site
#        if 'US-SRM' in df_param_sub.Site.unique():
#            param_SRM = df_param_sub[df_param_sub.Site == 'US-SRM']
#            param_SRM = param_SRM[~param_SRM.Param_ID.str.contains('__11')] # remove PFT=11 (15%) for US-SRM
#            df_param_sub = df_param_sub[~df_param_sub.Site.str.contains('SRM')]
#            df_param_sub = pd.concat([df_param_sub,param_SRM], ignore_index=True)
#        
#        # choose the parameters for PFT=05 (TeBE=55%) for the US-Ses site
#        if 'US-Ses' in df_param_sub.Site.unique():
#            param_Ses = df_param_sub[df_param_sub.Site == 'US-Ses']
#            param_Ses = param_Ses[~param_Ses.Param_ID.str.contains('__11')] # remove PFT=11 (25%) for US-Ses
#            df_param_sub = df_param_sub[~df_param_sub.Site.str.contains('Ses')]
#            df_param_sub = pd.concat([df_param_sub,param_Ses], ignore_index=True)
#        
#        # choose the parameters for PFT=04 (TeNE=60%) for the US-SRM site
#        if 'US-Mpj' in df_param_sub.Site.unique():
#            param_Mpj = df_param_sub[df_param_sub.Site == 'US-Mpj']
#            param_Mpj = param_Mpj[~param_Mpj.Param_ID.str.contains('__11')] # remove PFT=11 (20%) for US-Mpj
#            df_param_sub = df_param_sub[~df_param_sub.Site.str.contains('Mpj')]
#            df_param_sub = pd.concat([df_param_sub,param_Mpj], ignore_index=True)
#        
#        df_param_plot = df_param_sub[df_param_sub['Paramtype'].isin(['Param_post','Param_error_min','Param_error_max','Param_default','Param_min','Param_max'])]
##        df_param_plot = df_param.loc[df_param['Paramtype'] != 'Param_error'] 
#        g = sns.relplot(x='Site', y='Paramvalue', hue='Paramtype', style='Paramtype', 
#                size = 'Paramtype', col='Param_name', 
#                data=df_param_plot, palette=['purple','grey','orange','orange','purple','purple'], 
#                col_wrap=5, height=12, aspect=0.5,
#                markers=["o","X","P","P",'^','v'], sizes=(30,50),
#                facet_kws=dict(sharey=False),
#                size_order = ['Param_post','Param_error_min','Param_error_max','Param_default','Param_min','Param_max'])  
#
#        g.set(xlabel='', ylabel='')
#        plt.subplots_adjust(top=0.91,left=0.035,bottom=0.1,right=0.99)
#        g._legend.set_bbox_to_anchor([0.93,0.17])
#        # replace labels
#        new_labels = ['','Posterior parameters','Prior parameters','Parameter bounds','Parameter bounds','Parameter uncertainties','Parameter uncertainties']
#        for t, l in zip(g._legend.texts, new_labels): t.set_text(l)
#            
#        g.fig.suptitle('Parameters ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i])   
#        g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='Parameter values', size=12, rotation=90)
#            
#        for nx, ax in enumerate(g.axes.flat):
##            nx = 0
##            ax = g.axes.flat[nx]
#            plt.setp(ax.get_xticklabels(), rotation=90)
#                
#            mask1 = (df_param_sub["Paramtype"]=="Param_post") & (df_param_sub["Param_name"]==df_param_sub['Param_name'].unique()[nx])
#            df1 = df_param_sub.loc[mask1, ["Paramvalue"]]
#            mask2 = (df_param_sub["Paramtype"]=="Param_error") & (df_param_sub["Param_name"]==df_param_sub['Param_name'].unique()[nx])
#            df2 = df_param_sub.loc[mask2, ["Paramvalue"]]
#            ax.errorbar(df_param_sub.Site.unique(), df1.values, yerr=df2.values, fmt='.', lw=0.3, capsize=0, 
#                            capthick=0, color='grey')    
#            
#            mask3 = (df_param_sub["Paramtype"]=="Param_min") & (df_param_sub["Param_name"]==df_param_sub['Param_name'].unique()[nx])
#            df3 = df_param_sub.loc[mask3, ["Paramvalue"]]
#            mask4 = (df_param_sub["Paramtype"]=="Param_max") & (df_param_sub["Param_name"]==df_param_sub['Param_name'].unique()[nx])
#            df4 = df_param_sub.loc[mask4, ["Paramvalue"]]
#            ax.set(ylim=(min(df3.Paramvalue)-abs(min(df3.Paramvalue)*0.05), max(df4.Paramvalue)+abs(max(df4.Paramvalue)*0.05)))
#        
#        g.fig.savefig(outdir + 'Parameters/' + 'Parameters_optim_' + optim_name[i] + '_Csink.png')
#        plt.close()
        

        


    ##%% plotting the ratio between the posterior and prior RMSE of daily time series
    ## ---------------------
    ## ---------------------
    df_rmse_ratio = df_rmse.loc[df_rmse['Datatype'] == 'Ratio']
    df_rmse_ratio = df_rmse_ratio.drop('Datatype', 1)
    
    g = sns.catplot(x='Flux', y='RMSE', col='Site', data=df_rmse_ratio, kind='bar', order=['NEE','GPP','Reco'],
                     palette=['grey','blue','orange'], height=12, aspect=0.5, col_order = site_names)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.91,left=0.03)
    # entire figure title
    g.fig.suptitle('Ratio between the posterior and prior RMSE of daily time series fit ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i])   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
                    
    # save figure
    g.fig.savefig(outdir + 'RMSE/' + 'RMSE_ratio_daily_optim_' + optim_name[i] + '_Csink.png')
    plt.close()
    
    
    ##%% plotting the Fluxes monthly correlation coefficents
    ## ---------------------
    ## ---------------------
    df_month_corr = df_month.copy()
    df_month_corr['Month_year'] = pd.to_datetime(df_month_corr['Date']).dt.to_period('M')
    df_month_corr.drop(['Month', 'Date'], axis=1, inplace=True)
    
    flux_monthly_corr = df_month_corr.groupby(['Site','Flux'])[['Obs','Prior','Post']].corr()
    flux_monthly_corr.drop(['Prior', 'Post'], axis=1, inplace=True)
    flux_monthly_corr.drop('Obs', level=2, inplace=True)
    flux_monthly_corr.reset_index(level=0, inplace=True)
    flux_monthly_corr.reset_index(level=0, inplace=True)
    flux_monthly_corr.reset_index(level=0, inplace=True)
    flux_monthly_corr.columns = ['Simulation','Flux','Site','Monthly_corr']
        
    flux_monthly_corr = pd.merge(df_site_names,flux_monthly_corr,left_on='Site',right_on='Site',how='outer')
    nee_corr_prior = flux_monthly_corr[(flux_monthly_corr['Flux'] == 'NEE') & (flux_monthly_corr['Simulation'] == 'Prior')]
    nee_corr_post = flux_monthly_corr[(flux_monthly_corr['Flux'] == 'NEE') & (flux_monthly_corr['Simulation'] == 'Post')]
    gpp_corr_prior = flux_monthly_corr[(flux_monthly_corr['Flux'] == 'GPP') & (flux_monthly_corr['Simulation'] == 'Prior')]
    gpp_corr_post = flux_monthly_corr[(flux_monthly_corr['Flux'] == 'GPP') & (flux_monthly_corr['Simulation'] == 'Post')]
    reco_corr_prior = flux_monthly_corr[(flux_monthly_corr['Flux'] == 'Reco') & (flux_monthly_corr['Simulation'] == 'Prior')]
    reco_corr_post = flux_monthly_corr[(flux_monthly_corr['Flux'] == 'Reco') & (flux_monthly_corr['Simulation'] == 'Post')]
    
    #monthly fluxes R plot, seperated by fluxes (according to Natasha's suggestion)
    fig = plt.figure(figsize = (15,5))
    ind = np.arange(len(site_names)) + 0.15
    width = 0.35
    xtra_space = 0.05
    
    ax1 = fig.add_subplot(1,3,1)
    rect1 = ax1.bar(ind, nee_corr_prior.Monthly_corr, width, color='green') 
    rect2 = ax1.bar(ind + width, nee_corr_post.Monthly_corr, width, color='red') 
    ax2 = fig.add_subplot(1,3,2)
    rect3 = ax2.bar(ind + xtra_space, gpp_corr_prior.Monthly_corr, width, color='green')
    rect4 = ax2.bar(ind + width + xtra_space, gpp_corr_post.Monthly_corr, width, color='red')
    ax3 = fig.add_subplot(1,3,3)
    rect5 = ax3.bar(ind + xtra_space, reco_corr_prior.Monthly_corr, width, color='green')
    rect6 = ax3.bar(ind + width + xtra_space, reco_corr_post.Monthly_corr, width, color='red')
        
    ax1.legend((rect1,rect2), ('Prior','Post'),loc='best',ncol=1)
    ax1.set_title('NEE Flux')
    ax2.set_title('Flux monthly wrt observation ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i]
        + '\n' + 'GPP Flux')
    ax3.set_title('$R_{eco}$ Flux')
    ax1.set_xticks(ind+width)
    ax2.set_xticks(ind+width)
    ax3.set_xticks(ind+width)
    ax1.set_xticklabels(site_names, rotation=45)
    ax2.set_xticklabels(site_names, rotation=45)
    ax3.set_xticklabels(site_names, rotation=45)
    axisYmid = (max(max(nee_corr_prior.Monthly_corr), max(nee_corr_post.Monthly_corr)) + 
                min(min(nee_corr_prior.Monthly_corr), min(nee_corr_post.Monthly_corr)))/2
    ax1.text(x=-2, y=axisYmid, verticalalignment='center', s='Correlation coefficient (R)', size=12, rotation=90)
    
    plt.savefig(outdir + 'Flux_monthly/' + 'Flux_monthly_corrcoef_optim_' + optim_name[i] + '_Csink.png')
    plt.close()

        
    ##%% plotting the Fluxes annual anomaly correlation coefficents
    ## ---------------------
    ## ---------------------
    flux_corr = anomaly_site_all.groupby(['Site','Flux'])[['Obs','Prior','Post']].corr()
    flux_corr.drop(['Prior', 'Post'], axis=1, inplace=True)
    flux_corr.drop('Obs', level=2, inplace=True)
    flux_corr.reset_index(level=0, inplace=True)
    flux_corr.reset_index(level=0, inplace=True)
    flux_corr.reset_index(level=0, inplace=True)
    flux_corr.columns = ['Simulation','Flux','Site','Anomaly']    
    
    ## plotting Fluxes annual anomaly correlation coefficents according to C sink
    ## ---------------------
    ## ---------------------
    flux_corr = pd.merge(df_site_names,flux_corr,left_on='Site',right_on='Site',how='outer')
    
    nee_corr_prior = flux_corr[(flux_corr['Flux'] == 'NEE') & (flux_corr['Simulation'] == 'Prior')]
    nee_corr_post = flux_corr[(flux_corr['Flux'] == 'NEE') & (flux_corr['Simulation'] == 'Post')]
    gpp_corr_prior = flux_corr[(flux_corr['Flux'] == 'GPP') & (flux_corr['Simulation'] == 'Prior')]
    gpp_corr_post = flux_corr[(flux_corr['Flux'] == 'GPP') & (flux_corr['Simulation'] == 'Post')]
    reco_corr_prior = flux_corr[(flux_corr['Flux'] == 'Reco') & (flux_corr['Simulation'] == 'Prior')]
    reco_corr_post = flux_corr[(flux_corr['Flux'] == 'Reco') & (flux_corr['Simulation'] == 'Post')]
    
    # plots seperated by fluxes (according to Natasha's suggestion)
    fig = plt.figure(figsize = (15,5))
    ind = np.arange(len(site_names)) + 0.15
    width = 0.35
    xtra_space = 0.05
    
    ax1 = fig.add_subplot(1,3,1)
    rect1 = ax1.bar(ind, nee_corr_prior.Anomaly, width, color='green') 
    rect2 = ax1.bar(ind + width, nee_corr_post.Anomaly, width, color='red') 
    ax2 = fig.add_subplot(1,3,2)
    rect3 = ax2.bar(ind + xtra_space, gpp_corr_prior.Anomaly, width, color='green')
    rect4 = ax2.bar(ind + width + xtra_space, gpp_corr_post.Anomaly, width, color='red')
    ax3 = fig.add_subplot(1,3,3)
    rect5 = ax3.bar(ind + xtra_space, reco_corr_prior.Anomaly, width, color='green')
    rect6 = ax3.bar(ind + width + xtra_space, reco_corr_post.Anomaly, width, color='red')
        
    ax1.legend((rect1,rect2), ('Prior','Post'),loc='best',ncol=1)
    ax1.set_title('NEE Flux')
    ax2.set_title('Flux annual annomaly wrt observation ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i]
        + '\n' + 'GPP Flux')
    ax3.set_title('$R_{eco}$ Flux')
    ax1.set_xticks(ind+width)
    ax2.set_xticks(ind+width)
    ax3.set_xticks(ind+width)
    ax1.set_xticklabels(site_names, rotation=45)
    ax2.set_xticklabels(site_names, rotation=45)
    ax3.set_xticklabels(site_names, rotation=45)
    axisYmid = (max(max(nee_corr_prior.Anomaly), max(nee_corr_post.Anomaly)) + 
                min(min(nee_corr_prior.Anomaly), min(nee_corr_post.Anomaly)))/2
    ax1.text(x=-2, y=axisYmid, verticalalignment='center', s='Correlation coefficient (R)', size=12, rotation=90)
    
    plt.savefig(outdir + 'Flux_annual_anomaly/' + 'Flux_annual_anomaly_corrcoef_optim_' + optim_name[i] + '_Csink.png')
    plt.close()
    
    # appending the Correlation coefficient (R) data frame with all optimizations
    flux_corr['optim'] = i+1
    flux_corr['optim_name'] = [j.lower().replace(',', '').split('_iter') for j in optim_name][i][0]
    flux_corr_all = flux_corr_all.append(flux_corr)
    
           
    ##%% plotting the Fluxes annual anomaly slope from the linear least-squares regression
    ## ---------------------
    ## ---------------------
#    from scipy import stats
    flux_slope = flux_corr
    flux_slope.Anomaly = np.nan
    
    for nf, sf in enumerate(flux_names):     
#        nf = 0
        sf = flux_names[nf]
        anomaly_site_all_flux = anomaly_site_all[anomaly_site_all['Flux'] == sf]
        for ns, ss in enumerate(site_names): 
#            ns = 0
            ss = site_names[ns]
            anomaly_site_all_flux_site = anomaly_site_all_flux[anomaly_site_all_flux['Site'] == ss]
            
            slope_prior = linregress(anomaly_site_all_flux_site.Obs, anomaly_site_all_flux_site.Prior)
            slope_post = linregress(anomaly_site_all_flux_site.Obs, anomaly_site_all_flux_site.Post)
            
            flux_slope.Anomaly[(flux_slope['Flux'] == sf) & (flux_slope['Site'] == ss) & (flux_slope['Simulation'] == 'Prior')] = slope_prior.slope
            flux_slope.Anomaly[(flux_slope['Flux'] == sf) & (flux_slope['Site'] == ss) & (flux_slope['Simulation'] == 'Post')] = slope_post.slope
    
    nee_slope_prior = flux_slope[(flux_slope['Flux'] == 'NEE') & (flux_slope['Simulation'] == 'Prior')]
    nee_slope_post = flux_slope[(flux_slope['Flux'] == 'NEE') & (flux_slope['Simulation'] == 'Post')]
    gpp_slope_prior = flux_slope[(flux_slope['Flux'] == 'GPP') & (flux_slope['Simulation'] == 'Prior')]
    gpp_slope_post = flux_slope[(flux_slope['Flux'] == 'GPP') & (flux_slope['Simulation'] == 'Post')]
    reco_slope_prior = flux_slope[(flux_slope['Flux'] == 'Reco') & (flux_slope['Simulation'] == 'Prior')]
    reco_slope_post = flux_slope[(flux_slope['Flux'] == 'Reco') & (flux_slope['Simulation'] == 'Post')]
    
    fig = plt.figure(figsize = (15,5))
    ind = np.arange(len(site_names)) + 0.15
    width = 0.35
    xtra_space = 0.05
    
    ax1 = fig.add_subplot(1,3,1)
    rect1 = ax1.bar(ind, nee_slope_prior.Anomaly, width, color='green') 
    rect2 = ax1.bar(ind + width, nee_slope_post.Anomaly, width, color='red') 
    ax2 = fig.add_subplot(1,3,2)
    rect3 = ax2.bar(ind + xtra_space, gpp_slope_prior.Anomaly, width, color='green')
    rect4 = ax2.bar(ind + width + xtra_space, gpp_slope_post.Anomaly, width, color='red')
    ax3 = fig.add_subplot(1,3,3)
    rect5 = ax3.bar(ind + xtra_space, reco_slope_prior.Anomaly, width, color='green')
    rect6 = ax3.bar(ind + width + xtra_space, reco_slope_post.Anomaly, width, color='red')
        
    ax1.legend((rect1,rect2), ('Prior','Post'),loc='best',ncol=1)
    ax1.set_title('NEE Flux')
    ax2.set_title('Flux annual annomaly wrt observation ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i]
        + '\n' + 'GPP Flux')
    ax3.set_title('$R_{eco}$ Flux')
    ax1.set_xticks(ind+width)
    ax2.set_xticks(ind+width)
    ax3.set_xticks(ind+width)
    ax1.set_xticklabels(site_names, rotation=45)
    ax2.set_xticklabels(site_names, rotation=45)
    ax3.set_xticklabels(site_names, rotation=45)
    
    axisYmid = (max(max(nee_slope_prior.Anomaly), max(nee_slope_post.Anomaly)) + 
                min(min(nee_slope_prior.Anomaly), min(nee_slope_post.Anomaly)))/2
    ax1.text(x=-2, y=axisYmid, verticalalignment='center', s='Slope of linear least-squares regression', size=12, rotation=90)
    
    plt.savefig(outdir + 'Flux_annual_anomaly/' + 'Flux_annual_anomaly_slope_optim_' + optim_name[i] + '_Csink.png')
    plt.close()

    # appending the Slope of regression dataframe for all optimizations
    flux_slope['optim'] = i+1
    flux_slope['optim_name'] = [j.lower().replace(',', '').split('_iter') for j in optim_name][i][0]
    flux_slope_all = flux_slope_all.append(flux_slope)
    
    
    
    ##%% make a table of all correlation coefficients for the daily, monthly and annual time series
    ## ---------------------
    ## ---------------------
    # create similar dataframe for daily data as monthly dataframes
    df_flux_pivot = df_flux.set_index(['Date','Site','Flux','variable'])['value'].unstack().reset_index()
    
    flux_daily_corr = df_flux_pivot.groupby(['Site','Flux'])[['Obs','Prior','Post']].corr()
    flux_daily_corr.drop(['Prior', 'Post'], axis=1, inplace=True)
    flux_daily_corr.drop('Obs', level=2, inplace=True)
    flux_daily_corr.reset_index(level=0, inplace=True)
    flux_daily_corr.reset_index(level=0, inplace=True)
    flux_daily_corr.reset_index(level=0, inplace=True)
    flux_daily_corr.columns = ['Simulation','Flux','Site','Daily_corr']
        
    flux_daily_corr = pd.merge(df_site_names,flux_daily_corr,left_on='Site',right_on='Site',how='outer')
    
    
    # create similar dataframe for yearly data as monthly dataframes
#    df_annual_pivot = df_annual.set_index(['Date','Site','Flux','variable'])['value'].unstack().reset_index()
    
    flux_annual_corr = df_annual.groupby(['Site','Flux'])[['Obs','Prior','Post']].corr()
    flux_annual_corr.drop(['Prior', 'Post'], axis=1, inplace=True)
    flux_annual_corr.drop('Obs', level=2, inplace=True)
    flux_annual_corr.reset_index(level=0, inplace=True)
    flux_annual_corr.reset_index(level=0, inplace=True)
    flux_annual_corr.reset_index(level=0, inplace=True)
    flux_annual_corr.columns = ['Simulation','Flux','Site','Yearly_corr']
        
    flux_annual_corr = pd.merge(df_site_names,flux_annual_corr,left_on='Site',right_on='Site',how='outer')
    
    
    #merge all three ddataframes
#    flux_daily_corr = flux_daily_corr.rename(columns={'Anomaly': 'Daily_sum'})
#    flux_monthly_corr = flux_monthly_corr.rename(columns={'Anomaly': 'Monthly_sum'})
#    flux_annual_corr = flux_corr.rename(columns={'Anomaly': 'Annual_sum'})
    flux_corr_total = pd.merge(flux_daily_corr,flux_monthly_corr,left_on=['Site','Simulation','Flux'],right_on=['Site','Simulation','Flux'],how='outer')
    flux_corr_total = pd.merge(flux_corr_total,flux_annual_corr,left_on=['Site','Simulation','Flux'],right_on=['Site','Simulation','Flux'],how='outer')
    flux_corr_total = flux_corr_total.round(3) # round up the numbers to 3 decimal places
    
    #write the dataframe in outputs directory
    flux_corr_total.to_csv(outdir + 'Flux_corr/' + 'Flux_corr_optim_' + optim_name[i] + '.csv', encoding='utf-8', index=False)
    
    
    ##%% make a table of all slope for the daily, monthly and annual time series
    ## ---------------------
    ## ---------------------
    # create similar dataframe for daily slope data as daily correlation coefficients dataframes
    flux_daily_slope = flux_daily_corr[['Simulation','Flux','Site']]
    flux_daily_slope['Daily_slope'] = 0
    
    for nf, sf in enumerate(flux_names):     
#        nf = 0
        sf = flux_names[nf]
        df_flux_pivot_slope = df_flux_pivot[df_flux_pivot['Flux'] == sf]
        for ns, ss in enumerate(site_names): 
#            ns = 0
            ss = site_names[ns]
            df_flux_pivot_slope_site = df_flux_pivot_slope[df_flux_pivot_slope['Site'] == ss]
            
            slope_prior = linregress(df_flux_pivot_slope_site.Obs, df_flux_pivot_slope_site.Prior)
            slope_post = linregress(df_flux_pivot_slope_site.Obs, df_flux_pivot_slope_site.Post)
            
            flux_daily_slope.Daily_slope[(flux_daily_slope['Flux'] == sf) & (flux_daily_slope['Site'] == ss) & (flux_daily_slope['Simulation'] == 'Prior')] = slope_prior.slope
            flux_daily_slope.Daily_slope[(flux_daily_slope['Flux'] == sf) & (flux_daily_slope['Site'] == ss) & (flux_daily_slope['Simulation'] == 'Post')] = slope_post.slope
    flux_daily_slope = pd.merge(df_site_names,flux_daily_slope,left_on='Site',right_on='Site',how='outer')
    

    # create similar dataframe for monthly slope data
    flux_monthly_slope = flux_monthly_corr[['Simulation','Flux','Site']]
    flux_monthly_slope['Monthly_slope'] = 0
    
    for nf, sf in enumerate(flux_names):     
#        nf = 0
        sf = flux_names[nf]
        df_flux_pivot_slope = df_month_corr[df_month_corr['Flux'] == sf]
        for ns, ss in enumerate(site_names): 
#            ns = 0
            ss = site_names[ns]
            df_flux_pivot_slope_site = df_flux_pivot_slope[df_flux_pivot_slope['Site'] == ss]
            
            slope_prior = linregress(df_flux_pivot_slope_site.Obs, df_flux_pivot_slope_site.Prior)
            slope_post = linregress(df_flux_pivot_slope_site.Obs, df_flux_pivot_slope_site.Post)
            
            flux_monthly_slope.Monthly_slope[(flux_monthly_slope['Flux'] == sf) & (flux_monthly_slope['Site'] == ss) & (flux_monthly_slope['Simulation'] == 'Prior')] = slope_prior.slope
            flux_monthly_slope.Monthly_slope[(flux_monthly_slope['Flux'] == sf) & (flux_monthly_slope['Site'] == ss) & (flux_monthly_slope['Simulation'] == 'Post')] = slope_post.slope
    flux_monthly_slope = pd.merge(df_site_names,flux_monthly_slope,left_on='Site',right_on='Site',how='outer')
    
    
    
    # create similar dataframe for yearly slope data
    flux_annual_slope = flux_annual_corr[['Simulation','Flux','Site']]
    flux_annual_slope['Yearly_slope'] = 0
    
    for nf, sf in enumerate(flux_names):     
#        nf = 0
        sf = flux_names[nf]
        df_flux_pivot_slope = df_annual[df_annual['Flux'] == sf]
        for ns, ss in enumerate(site_names): 
#            ns = 0
            ss = site_names[ns]
            df_flux_pivot_slope_site = df_flux_pivot_slope[df_flux_pivot_slope['Site'] == ss]
            
            slope_prior = linregress(df_flux_pivot_slope_site.Obs, df_flux_pivot_slope_site.Prior)
            slope_post = linregress(df_flux_pivot_slope_site.Obs, df_flux_pivot_slope_site.Post)
            
            flux_annual_slope.Yearly_slope[(flux_annual_slope['Flux'] == sf) & (flux_annual_slope['Site'] == ss) & (flux_annual_slope['Simulation'] == 'Prior')] = slope_prior.slope
            flux_annual_slope.Yearly_slope[(flux_annual_slope['Flux'] == sf) & (flux_annual_slope['Site'] == ss) & (flux_annual_slope['Simulation'] == 'Post')] = slope_post.slope
    flux_annual_slope = pd.merge(df_site_names,flux_annual_slope,left_on='Site',right_on='Site',how='outer')
    
    
    #merge all three ddataframes
#    flux_daily_corr = flux_daily_corr.rename(columns={'Anomaly': 'Daily_sum'})
#    flux_monthly_corr = flux_monthly_corr.rename(columns={'Anomaly': 'Monthly_sum'})
#    flux_annual_corr = flux_corr.rename(columns={'Anomaly': 'Annual_sum'})
    flux_slope_total = pd.merge(flux_daily_slope,flux_monthly_slope,left_on=['Site','Simulation','Flux'],right_on=['Site','Simulation','Flux'],how='outer')
    flux_slope_total = pd.merge(flux_slope_total,flux_annual_slope,left_on=['Site','Simulation','Flux'],right_on=['Site','Simulation','Flux'],how='outer')
    flux_slope_total = flux_slope_total.round(3) # round up the numbers to 3 decimal places
    
    #write the dataframe in outputs directory
    flux_slope_total.to_csv(outdir + 'Flux_corr/' + 'Flux_slope_optim_' + optim_name[i] + '.csv', encoding='utf-8', index=False)
    
    
    
    
    
    ##%% plotting the MSD
    ## ---------------------
    ## ---------------------
    # plotting the MSD and MSD reduction with seperated sites (Natasha's suggestion)
    ind = 0.5
    width = 0.4
    for nf, sf in enumerate(flux_names):     
#        nf = 0
        sf = flux_names[nf]
        df_msd_flux = df_msd[df_msd['Flux'] == sf]
        df_msd_flux.drop('optim', axis=1, inplace=True)
        df_msd_flux.drop('optim_name', axis=1, inplace=True)
        
        fig, ax = plt.subplots(nrows=2, ncols=len(site_names), figsize=(10,10))
        for ns, ss in enumerate(site_names): 
#            ns = 0
            ss = site_names[ns]
            df_msd_site = df_msd_flux[df_msd_flux['Site'] == ss]
            df_msd_site.drop(['mse','r','Flux'], axis=1, inplace=True)
            df_msd_site.set_index('Site', inplace=True)
            df_msd_site_reduction = df_msd_site.copy()
            df_msd_site_reduction.drop('Sim', axis=1, inplace=True)
            df_msd_site_reduction = df_msd_site_reduction.iloc[:-1]
            
            df_msd_site_reduction.set_value(ss, 'bias', 1 - (df_msd_site.bias[df_msd_site['Sim'] == 'post'] / df_msd_site.bias[df_msd_site['Sim'] == 'prior']))
            df_msd_site_reduction.set_value(ss, 'variance', 1 - (df_msd_site.variance[df_msd_site['Sim'] == 'post'] / df_msd_site.variance[df_msd_site['Sim'] == 'prior']))
            df_msd_site_reduction.set_value(ss, 'phase', 1 - (df_msd_site.phase[df_msd_site['Sim'] == 'post'] / df_msd_site.phase[df_msd_site['Sim'] == 'prior']))
            
            df_msd_site_post = df_msd_site[df_msd_site['Sim'] == 'post']
            df_msd_site_prior = df_msd_site[df_msd_site['Sim'] == 'prior']
            df_msd_site_post.variance = df_msd_site_post.variance + df_msd_site_post.bias
            df_msd_site_post.phase = df_msd_site_post.variance + df_msd_site_post.phase
            df_msd_site_prior.variance = df_msd_site_prior.variance + df_msd_site_prior.bias
            df_msd_site_prior.phase = df_msd_site_prior.variance + df_msd_site_prior.phase
            
            rect1 = ax[0][ns].bar(ind, df_msd_site_prior.phase, width, color='bisque', hatch='.') 
            rect2 = ax[0][ns].bar(ind, df_msd_site_prior.variance, width, color='lime', hatch='.') 
            rect3 = ax[0][ns].bar(ind, df_msd_site_prior.bias, width, color='mediumorchid', hatch='.')
            
            xtra_space = 0.05
            rect4 = ax[0][ns].bar(ind + width + xtra_space, df_msd_site_post.phase, width, color='orange', hatch='/') 
            rect5 = ax[0][ns].bar(ind + width + xtra_space, df_msd_site_post.variance, width, color='green', hatch='/') 
            rect6 = ax[0][ns].bar(ind + width + xtra_space, df_msd_site_post.bias, width, color='purple', hatch='/')
            
            
            rect7 = ax[1][ns].bar(np.arange(3), df_msd_site_reduction.values[0], color=['purple','green','orange']) 
            
            if ns == 0:
                axisY1mid = max(df_msd_site_prior.phase.values,df_msd_site_post.phase.values)/2
                if max(df_msd_site_reduction.values[0]) < 0:
                    axisY2mid = min(df_msd_site_reduction.values[0])/2
                elif min(df_msd_site_reduction.values[0]) < 0 and max(df_msd_site_reduction.values[0]) > 0:
                    axisY2mid = (max(df_msd_site_reduction.values[0]) + min(df_msd_site_reduction.values[0]))/2
                elif max(df_msd_site_reduction.values[0]) > 0:
                    axisY2mid = max(df_msd_site_reduction.values[0])/2
                
            ax[0][ns].set_xticks([])
            ax[1][ns].set_xticks([])
            ax[1][ns].set_xlabel(ss)
        
            ax[1][ns].axhline(y=0, linestyle='-', color='grey') # add a horizontal line at y=0
        
        # add lagends
        ax[0][len(site_names)-1].legend((rect3,rect2,rect1,rect6,rect5,rect4), ('Prior $Bias^{2}$','Prior $Variance^{2}$','Prior Phase','Post $Bias^{2}$','Post $Variance^{2}$','Post Phase'), 
          bbox_to_anchor=(1.2, 1.25), ncol=2)
        ax[1][len(site_names)-1].legend((rect7), ('$Bias^{2}$','$Variance^{2}$','Phase'), bbox_to_anchor=(1.2, 1.12), ncol=3)
    
        # add labels, title and axes ticks
        ax[0][round((len(site_names)-1)/2)].set_title(fig_title[nf]+' Flux')
            
        plt.subplots_adjust(wspace=0.8)
    
        #add figure title
        fig.suptitle('Mean square deviation (MSD) decomposition ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i]) 
           
        ax[0][0].text(x=-0.7, y=axisY1mid, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase ($(gCm^{-2}d^{-1})^{2}$)', size=12, rotation=90)
        ax[1][0].text(x=-3.7, y=axisY2mid, verticalalignment='center', s='Reduction in $Bias^{2}$ / $Variance^{2}$ / Phase', size=12, rotation=90)
            
        plt.savefig(outdir + 'MSD/' + 'Mean_square_deviation_decomposition_' + flux_names[nf] + '_optim_' + optim_name[i] + '_Csink.png')
        plt.close()
    
    
    ##%% plotting the mean monthly flux time series for all sites
    ## ---------------------
    ## ---------------------
    #get the mean monthly fluxes cycle
    g = sns.FacetGrid(df_month_mean, col="Site", hue="Datatype", row="Flux", margin_titles=True,
                       col_order = site_names, palette=pal, sharey=False, height=12, aspect=0.5)
    g = g.map(plt.plot, "Month", "Monthly_mean", alpha=.7)
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.91,left=0.035)
    plt.legend(bbox_to_anchor=(1, 3.45), ncol=3)
    g.set(xticks=df_month_mean.Month.unique()[1::2])
    
    g.fig.suptitle('Mean monthly fluxes ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i])   
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='NEE / GPP / $R_{eco}$ ($gCm^{-2}month^{-1}$)', size=12, rotation=90)
    g.fig.text(x=0.5, y=0.02, verticalalignment='center', s='Months', size=12, rotation=0)
    
    #draw a horizontal line at NEE = 0 to see the C sourse / sink
    for nx, ax in enumerate(g.axes.flat):
        if nx <= len(site_names):
            ax.axhline(y=0, linestyle='--', color='grey') # add a horizontal line at NEE = 0
    
    g.fig.savefig(outdir + 'Mean_monthly_fluxes/' + 'Mean_monthly_fluxes_optim_' + optim_name[i] + '_Csink.png')
    plt.close()

    
    ##%% plotting the flux time series for all sites
    ## ---------------------
    ## ---------------------
    # extract all 3 flux data into dataframe for each site
    plot_flux_text = [plot_nee_text,plot_gpp_text,plot_reco_text]
    
    for nf, sf in enumerate(flux_names):
#        nf = 0
        sf = flux_names[nf]
        
        df_ind_flux = df_flux.loc[df_flux['Flux'] == sf]
        
        sns.set(font_scale=1.2, style='white')
        g = sns.FacetGrid(df_ind_flux, col="Site", hue="variable", palette=pal, col_wrap=col_wrap, sharex=False, 
                      sharey=False, height=12, aspect=0.5)
        g = g.map(plt.plot, "Date", "value", alpha=.7)
        plt.legend(loc='best', ncol=3)
        g.set(xlabel='', ylabel='')
        plt.subplots_adjust(top=0.91,left=0.04,bottom=0.035,right=0.995)
        
        # add title with RMSE values      
        for ax, title in zip(g.axes.flat,plot_flux_text[nf]):
            ax.set_title(title)
        # entire figure title
#        g.fig.suptitle(fig_title[nf] + ' observation and model simulations ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i])   
        # one ylabel
        g.fig.text(x=0.005, y=0.5, verticalalignment='center', s=sf+' ($gCm^{-2}day^{-1}$)', size=12, rotation=90)
        
        # save figure
        g.fig.savefig(outdir + 'Flux_timeseries/' + sf + '_optim_' + optim_name[i] + '_Csink.png')
        plt.close()
        
        
        ## ---------------------
        ## ---------------------
        # Draw a barplot to show the annual flux anomaly
        anomaly_ind_flux = df_anomaly.loc[df_anomaly['Flux'] == sf]
        
        g = sns.catplot(x='Year', y='Anomaly', hue='Datatype', 
                       col='Site', data=anomaly_ind_flux, kind='bar', col_wrap=col_wrap, col_order = site_names,
                       sharex=False, sharey=False, palette=pal, height=12, aspect=0.5, legend=False)    
        g.set(xlabel='', ylabel='')
        plt.legend(loc='best', ncol=3)
        plt.subplots_adjust(top=0.91,left=0.04)
        # entire figure title
        g.fig.suptitle(fig_title[nf] + ' annomaly ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i])   
        # one ylabel
        g.fig.text(x=0.005, y=0.5, verticalalignment='center', s=fig_title[nf]+' anomaly ($gCm^{-2}year^{-1}$)', size=12, rotation=90)
                
        # add title with RMSE values      
        for ax, title in zip(g.axes.flat,plot_flux_text[nf]):
            ax.set_title(title)
        
        # save figure
        g.fig.savefig(outdir + 'Flux_annual_anomaly/' + sf + '_annomaly_optim_' + optim_name[i] + '_Csink.png')
        plt.close()
    

#%% 
# - load flux corr data for Opt1 P1
flux_corr_nee_p1 = pd.read_csv(outdir + 'Flux_corr/' + 'Flux_corr_optim_' + optim_name[7] + '.csv', sep=",")
    
flux_corr_nee_p1_melt = pd.melt(flux_corr_nee_p1, id_vars=('Site','Simulation','Flux'), var_name="Frequency", value_name="Corr")            

sns.set(font_scale=1, style='white')   
for nf, sf in enumerate(flux_names):
#        nf = 0
    sf = flux_names[nf]
        
    #Box-whisker plots for correlation coefficient
    flux_corr_plot = flux_corr_nee_p1_melt.loc[flux_corr_nee_p1_melt['Flux'] == sf]
#        flux_corr_plot.optim[flux_corr_plot.Simulation == 'Prior'] = 0
#        flux_corr_plot.optim_name[flux_corr_plot.Simulation == 'Prior'] = 'Prior'
        
    g = sns.boxplot(x='Frequency', y='Corr', hue='Simulation', data=flux_corr_plot, dodge=True)  
    g.figure.set_size_inches(3,4)
    g.set(xlabel='', ylabel='')
        
    if nf==0: 
        plt.subplots_adjust(top=0.93,bottom=0.08,left=0.22,right=0.99)
        plt.text(x=-1.3, y=(max(flux_corr_plot.Corr)+min(flux_corr_plot.Corr))/2, verticalalignment='center', s='Correlation coefficient (R)', size=14, rotation=90)
    else: plt.subplots_adjust(top=0.93,bottom=0.08,left=0.15,right=0.99)
    plt.legend(title='Optimization', bbox_to_anchor=(0.55, 0.2), ncol=1, fontsize=8)
    plt.suptitle(subplot_title[nf]+' '+sf)   
    plt.ylim(-0.5, 1)
    
    g.set_xticklabels(['Daily','Monthly','Yearly'], fontsize=10)
    new_labels = ['Prior','P1']
    leg = g.axes.get_legend()
    for t, l in zip(leg.texts, new_labels): t.set_text(l)
        
    if nf > 0: 
        g.legend_.remove()
        
    plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_P1_anomaly_corrcoef_box_whisker.png')
    plt.close()

#%% 
# - load flux corr data for Opt1 P1
flux_slope_nee_p1 = pd.read_csv(outdir + 'Flux_corr/' + 'Flux_slope_optim_' + optim_name[7] + '.csv', sep=",")
    
flux_slope_nee_p1_melt = pd.melt(flux_slope_nee_p1, id_vars=('Site','Simulation','Flux'), var_name="Frequency", value_name="Slope")            

sns.set(font_scale=1, style='white')   
for nf, sf in enumerate(flux_names):
#        nf = 1
    sf = flux_names[nf]
        
    #Box-whisker plots for correlation coefficient
    flux_slope_plot = flux_slope_nee_p1_melt.loc[flux_slope_nee_p1_melt['Flux'] == sf]        
    g = sns.boxplot(x='Frequency', y='Slope', hue='Simulation', data=flux_slope_plot, dodge=True)  
    g.figure.set_size_inches(3,4)
    g.set(xlabel='', ylabel='')
    g.legend_.remove()
    
    if nf==0: 
        plt.subplots_adjust(top=0.93,bottom=0.08,left=0.22,right=0.99)
        plt.text(x=-1.3, y=(max(flux_slope_plot.Slope)+min(flux_slope_plot.Slope))/2, verticalalignment='center', s='Slope of linear least-square regression', size=14, rotation=90)
    else: plt.subplots_adjust(top=0.93,bottom=0.08,left=0.17,right=0.99)
    plt.suptitle(subplot_title[nf+3]+' '+sf)       
    plt.ylim(-0.5, 1.5)
    g.set_xticklabels(['Daily','Monthly','Yearly'], fontsize=10)
        
    plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_P1_anomaly_slope_box_whisker.png')
    plt.close()


#%% plotting all mean monthly seasonal time series for all GPP+Reco optimizations
## ---------------------
## ---------------------
df_month_mean_all_plot = df_month_mean_all[df_month_mean_all['optim_name'].isin(['gpp_reco'])]
df_month_mean_all_plot.optim[df_month_mean_all_plot.Datatype == 'Obs'] = 'Obs'
df_month_mean_all_plot.optim[df_month_mean_all_plot.Datatype == 'Prior'] = 'Prior'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 1] = 'P1'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 2] = 'P2'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 3] = 'P3'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 4] = 'P4'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 5] = 'P5'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 6] = 'P6'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 7] = 'P7'
    
df_month_mean_all_plot = df_month_mean_all_plot.drop_duplicates(subset=['Month','Datatype','Site','Flux','optim'])
         
#pal1 = sns.cubehelix_palette(10, start=.5, rot=-.75)
#pal1[1] = "red"
#pal1[0] = "black"
#pal1 = ['black','red','magenta','royalblue','deepskyblue','blue','green','lime','cyan']
pal1 = ['black','red','magenta','royalblue','blue','green','lime','cyan','gold']
    
sns.set(font_scale=1.1, style='white')
g = sns.FacetGrid(df_month_mean_all_plot, col="Site", row='Flux', hue="optim", margin_titles=True,
                           col_order = site_names, palette=pal1, sharey=False, size=3, aspect=0.6)
g = g.map(plt.plot, "Month", "Monthly_mean", alpha=.7)
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.set(xlabel='', ylabel='')
plt.subplots_adjust(top=0.88,left=0.045)
plt.legend(bbox_to_anchor=(1, 3.7), ncol=3)
g.set(xticks=df_month_mean.Month.unique()[1::2])
plt.subplots_adjust(hspace=0.1, wspace=0.3)
       
g.fig.suptitle('Mean monthly fluxes ' + r"$\bf{FOR}$" + ' GPP+Reco optimizations')   
g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='NEE / GPP / $R_{eco}$ ($gCm^{-2}month^{-1}$)', size=12, rotation=90)
g.fig.text(x=0.5, y=0.02, verticalalignment='center', s='Months', size=12, rotation=0)
        
#draw a horizontal line at NEE = 0 to see the C sourse / sink
for nx, ax in enumerate(g.axes.flat):
    if nx < len(site_names):
        ax.axhline(y=0, linestyle='--', color='grey') # add a horizontal line at NEE = 0
    
g.fig.savefig(outdir + 'Mean_monthly_fluxes/' + 'Mean_monthly_fluxes_optim_GPP+Reco_all.png')
plt.close()


#%% plotting all mean monthly seasonal time series for all NEE optimizations
## ---------------------
## ---------------------
df_month_mean_all_plot = df_month_mean_all[df_month_mean_all['optim_name'].isin(['nee'])]
df_month_mean_all_plot.optim[df_month_mean_all_plot.Datatype == 'Obs'] = 'Obs'
df_month_mean_all_plot.optim[df_month_mean_all_plot.Datatype == 'Prior'] = 'Prior'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 8] = 'P1'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 9] = 'P2'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 10] = 'P3'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 11] = 'P4'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 12] = 'P5'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 13] = 'P6'
df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 14] = 'P7'
    
df_month_mean_all_plot = df_month_mean_all_plot.drop_duplicates(subset=['Month','Datatype','Site','Flux','optim'])
         
pal1 = ['black','red','blue','royalblue','green','lime','cyan','gold','magenta']
    
sns.set(font_scale=1.1, style='white')
g = sns.FacetGrid(df_month_mean_all_plot, col="Site", row='Flux', hue="optim", margin_titles=True,
                           col_order = site_names, palette=pal1, sharey=False, 
                           size=3, aspect=0.6)
g = g.map(plt.plot, "Month", "Monthly_mean", alpha=.7)
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.set(xlabel='', ylabel='')
plt.subplots_adjust(top=0.88,left=0.045)
plt.legend(bbox_to_anchor=(1, 3.7), ncol=3)
g.set(xticks=df_month_mean.Month.unique()[1::2])
        
g.fig.suptitle('Mean monthly fluxes ' + r"$\bf{FOR}$" + ' NEE optimizations')   
g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='NEE / GPP / $R_{eco}$ ($gCm^{-2}month^{-1}$)', size=12, rotation=90)
g.fig.text(x=0.5, y=0.02, verticalalignment='center', s='Months', size=12, rotation=0)
plt.subplots_adjust(hspace=0.1, wspace=0.3)
       
#draw a horizontal line at NEE = 0 to see the C sourse / sink
for nx, ax in enumerate(g.axes.flat):
    if nx < len(site_names):
        ax.axhline(y=0, linestyle='--', color='grey') # add a horizontal line at NEE = 0
    
g.fig.savefig(outdir + 'Mean_monthly_fluxes/' + 'Mean_monthly_fluxes_optim_NEE_all.png')
plt.close()


#plotting all mean monthly seasonal time series for all NEE optimizations with site type
#assign site type (source / pivot / sink)
df_month_mean_all_plot['site_type'] = 'Sink'
#    flux_msd_all.loc[flux_msd_all['Site'] == 'US-SO2', 'site_type'] = 'Pivot'
df_month_mean_all_plot.loc[df_month_mean_all_plot['Site'] == 'US-Wkg', 'site_type'] = 'Pivot'
df_month_mean_all_plot.loc[df_month_mean_all_plot['Site'] == 'US-Whs', 'site_type'] = 'Pivot'
df_month_mean_all_plot.loc[df_month_mean_all_plot['Site'] == 'US-SRG', 'site_type'] = 'Pivot'
df_month_mean_all_plot.loc[df_month_mean_all_plot['Site'] == 'US-SRM', 'site_type'] = 'Pivot'
df_month_mean_all_plot.loc[df_month_mean_all_plot['Site'] == 'US-Seg', 'site_type'] = 'Pivot'
    
#    flux_msd_all.loc[flux_msd_all['Site'] == 'US-Seg', 'site_type'] = 'Source'
df_month_mean_all_plot.loc[df_month_mean_all_plot['Site'] == 'US-Aud', 'site_type'] = 'Source'

#mean across site type
df_month_plot_mean = df_month_mean_all_plot.groupby(['Month','Datatype','Flux','optim','optim_name','site_type'], as_index=False).mean()

sns.set(font_scale=1.2, style='white')
g = sns.FacetGrid(df_month_plot_mean, row="site_type", col='Flux', hue="optim", margin_titles=True,
                           row_order = ['Sink','Pivot','Source'], col_order = ['NEE','GPP','Reco'], 
                           hue_order = ['Obs','Prior','P1','P2','P3','P4','P5','P6','P7'], palette=pal1, sharey=False, 
                           size=3, aspect=1.5)
g = g.map(plt.plot, "Month", "Monthly_mean", alpha=.7)
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.set(xlabel='', ylabel='')
plt.subplots_adjust(top=0.87,left=0.065)
plt.legend(bbox_to_anchor=(1, 3.75), ncol=3)
g.set(xticks=df_month_mean.Month.unique()[1::2])
        
#g.fig.suptitle('Mean monthly fluxes ' + r"$\bf{FOR}$" + ' NEE optimizations')   
g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='NEE / GPP / $R_{eco}$ ($gCm^{-2}month^{-1}$)', size=16, rotation=90)
g.fig.text(x=0.5, y=0.02, verticalalignment='center', s='Months', size=16, rotation=0)
plt.subplots_adjust(hspace=0.1, wspace=0.12)
       
#draw a horizontal line at NEE = 0 to see the C sourse / sink
for nx, ax in enumerate(g.axes.flat):
    if nx in {0,3,6}:
        ax.axhline(y=0, linestyle='--', color='grey') # add a horizontal line at NEE = 0

g.fig.text(x=0.14, y=0.85, verticalalignment='center', s='(a)', size=16)        
g.fig.text(x=0.45, y=0.85, verticalalignment='center', s='(b)', size=16)        
g.fig.text(x=0.75, y=0.85, verticalalignment='center', s='(c)', size=16)        
g.fig.text(x=0.14, y=0.58, verticalalignment='center', s='(d)', size=16)        
g.fig.text(x=0.45, y=0.58, verticalalignment='center', s='(e)', size=16)        
g.fig.text(x=0.75, y=0.58, verticalalignment='center', s='(f)', size=16)        
g.fig.text(x=0.14, y=0.31, verticalalignment='center', s='(g)', size=16)        
g.fig.text(x=0.45, y=0.31, verticalalignment='center', s='(h)', size=16)        
g.fig.text(x=0.75, y=0.31, verticalalignment='center', s='(i)', size=16)        
            
g.fig.savefig(outdir + 'Mean_monthly_fluxes/' + 'Mean_monthly_fluxes_optim_NEE_sitetype.png')
plt.close()

    
## plotting all mean monthly seasonal time series for P1 NEE optimizations
df_month_mean_all_plot_P1 = df_month_mean_all_plot[df_month_mean_all_plot['optim'].isin(['Obs','Prior','P1'])]
df_month_mean_all_plot_P1.optim[df_month_mean_all_plot_P1.optim == 'P1'] = 'Posterior'

sns.set(font_scale=1, style='white')
g = sns.FacetGrid(df_month_mean_all_plot_P1, col="Site", row='Flux', hue="optim", margin_titles=True,
                           col_order = site_names, palette=['black','red','blue'], sharey=False, 
                           size=3, aspect=0.6)
g = g.map(plt.plot, "Month", "Monthly_mean", alpha=.7)
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.set(xlabel='', ylabel='')
plt.subplots_adjust(top=0.92,left=0.045)
plt.legend(bbox_to_anchor=(1, 3.5), ncol=3)
g.set(xticks=df_month_mean.Month.unique()[1::2])
        
#g.fig.suptitle('Mean monthly fluxes ' + r"$\bf{FOR}$" + ' NEE optimizations')   
g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='NEE / GPP / $R_{eco}$ ($gCm^{-2}month^{-1}$)', size=12, rotation=90)
g.fig.text(x=0.5, y=0.02, verticalalignment='center', s='Months', size=12, rotation=0)
plt.subplots_adjust(hspace=0.1, wspace=0.35)
       
#draw a horizontal line at NEE = 0 to see the C sourse / sink
for nx, ax in enumerate(g.axes.flat):
    if nx < len(site_names):
        ax.axhline(y=0, linestyle='--', color='grey') # add a horizontal line at NEE = 0
    
g.fig.savefig(outdir + 'Mean_monthly_fluxes/' + 'Mean_monthly_fluxes_optim_NEE_P1.png')
plt.close()   
 

   
for nf, sf in enumerate(flux_names):
#    nf = 0
    sf = flux_names[nf]
    df_month_mean_all_plot = df_month_mean_all.loc[df_month_mean_all['Flux'] == sf]
    
    df_month_mean_all_plot = df_month_mean_all_plot[df_month_mean_all_plot['optim_name'].isin(['nee'])]
    df_month_mean_all_plot.optim[df_month_mean_all_plot.Datatype == 'Obs'] = 'Obs'
    df_month_mean_all_plot.optim[df_month_mean_all_plot.Datatype == 'Prior'] = 'Prior'
    df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 8] = 'P1'
    df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 9] = 'P2'
    df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 10] = 'P3'
    df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 11] = 'P4'
    df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 12] = 'P5'
    df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 13] = 'P6'
    df_month_mean_all_plot.optim[df_month_mean_all_plot.optim == 14] = 'P7'
    
    df_month_mean_all_plot = df_month_mean_all_plot.drop_duplicates(subset=['Month','Datatype','Site','Flux','optim'])
         
    sns.set(style='white')
    g = sns.FacetGrid(df_month_mean_all_plot, col="Site", hue="optim", margin_titles=True,
                           col_order = site_names, col_wrap=4, palette=pal1, sharey=False, 
                           size=2.5, aspect=1)
    g = g.map(plt.plot, "Month", "Monthly_mean", alpha=.7)
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.89,left=0.07)
    plt.legend(bbox_to_anchor=(1.1, 4), ncol=5)
    g.set(xticks=df_month_mean.Month.unique()[1::2])
        
#    g.fig.suptitle('Mean monthly ' + sf + r" $\bf{FOR}$ " + 'NEE optimization')   
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s=fig_title[nf] + ' ($gCm^{-2}month^{-1}$)', size=12, rotation=90)
#    g.fig.suptitle('Mean monthly fluxes ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i])   
#    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='NEE / GPP / $R_{eco}$ ($gCm^{-2}month^{-1}$)', size=12, rotation=90)
    g.fig.text(x=0.5, y=0.02, verticalalignment='center', s='Months', size=12, rotation=0)
        
    #draw a horizontal line at NEE = 0 to see the C sourse / sink
    if nf == 0:
        for nx, ax in enumerate(g.axes.flat):
#            if nx < len(site_names):
            ax.axhline(y=0, linestyle='--', color='grey') # add a horizontal line at NEE = 0
    
    g.fig.savefig(outdir + 'Mean_monthly_fluxes/' + sf + '_mean_monthly_fluxes_optim_NEE_all.png')
    plt.close()


    ## plotting all mean monthly seasonal time series for P1 NEE optimizations
    df_month_mean_P1 = df_month_mean_all_plot[df_month_mean_all_plot['optim'].isin(['Obs','Prior','P1'])]
    df_month_mean_P1.optim[df_month_mean_P1.optim == 'P1'] = 'Posterior'
    
    sns.set(style='white')
    g = sns.FacetGrid(df_month_mean_P1, col="Site", hue="optim", margin_titles=True,
                           col_order = site_names, col_wrap=4, palette=['black','red','blue'], sharey=False, 
                           size=2.5, aspect=1)
    g = g.map(plt.plot, "Month", "Monthly_mean", alpha=.7)
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.93,left=0.07)
    plt.legend(bbox_to_anchor=(1.1, 3.85), ncol=3)
    g.set(xticks=df_month_mean.Month.unique()[1::2])
        
#    g.fig.suptitle('Mean monthly ' + sf + r" $\bf{FOR}$ " + 'NEE optimization')   
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s=fig_title[nf] + ' ($gCm^{-2}month^{-1}$)', size=12, rotation=90)
#    g.fig.suptitle('Mean monthly fluxes ' + r"$\bf{FOR}$ " + optim[i] + ' optimization ' + r"$\bf{WITH}$ " + param_title[i])   
#    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='NEE / GPP / $R_{eco}$ ($gCm^{-2}month^{-1}$)', size=12, rotation=90)
    g.fig.text(x=0.5, y=0.02, verticalalignment='center', s='Months', size=12, rotation=0)
        
    #draw a horizontal line at NEE = 0 to see the C sourse / sink
    if nf == 0:
        for nx, ax in enumerate(g.axes.flat):
#            if nx < len(site_names):
            ax.axhline(y=0, linestyle='--', color='grey') # add a horizontal line at NEE = 0
    
    g.fig.savefig(outdir + 'Mean_monthly_fluxes/' + sf + '_mean_monthly_fluxes_optim_NEE_P1.png')
    plt.close()


#%% plotting all parameter uncertainty reductions and deviations
## ---------------------
## ---------------------
if param_plot == True: 
    #pd.set_option('display.max_rows', 50) # display upto 500 rows
    #pd.set_option('display.max_columns', 1000) # display upto 500 columns
    ## plotting all parameter uncertainty reductions
    ## ---------------------
    ## ---------------------
    boxplot_order=[8,9,10,11,12,13,14,1,2,3,4,5,6,7,15,16,17,18,19,20]
    df_param_group = pd.DataFrame(param_group)
    df_param_group['opt_param'] = df_param_group["opt"] + df_param_group["param"]

    df_param_plot = df_param_all.loc[df_param_all['Paramtype'] == 'uncer_reduc']
    df_param_plot = df_param_plot.drop('Paramtype', 1)
    df_param_plot = df_param_plot.drop('Param_ID', 1)
    df_param_plot = df_param_plot.drop('Site', 1)
    
    df_param_plot['uncer_reduc'] = df_param_plot.groupby(['Param_name','optim']).Paramvalue.transform('median')
    df_param_plot = df_param_plot.drop_duplicates(subset=['Param_name','optim'])
    
    # read param list to add the param types
    param_list = pd.read_excel(param_file)
    
    # merge the data frames to have param type
    df_param_plot = pd.merge(df_param_plot, param_list, how='right', on='Param_name')
#    df_param_plot = pd.merge(df_param_plot, param_list, left_on='Param_name',right_on='Param_name',how='outer')
    df_param_plot = df_param_plot.sort_values(['Param_type','Param_name'], ascending=True)
    df_param_plot.uncer_reduc = 100*df_param_plot.uncer_reduc
    
    df_param_xt = df_param_plot.copy() # mark the parameters for NEE opt with all param
    df_param_xt = df_param_xt.drop('Paramvalue', 1)
    
    #plot parameter deviation bar plot
    sns.set(font_scale=1, style='white')
    if any('paramP7' in s for s in optim_name): 
        g = sns.catplot(x='Param_name', y='uncer_reduc', row='optim', hue='Param_type', data=df_param_plot,
                    row_order=boxplot_order, kind='bar', size=3, aspect=4, margin_titles=True, legend=False, dodge=False) 
    else:
        g = sns.catplot(x='Param_name', y='uncer_reduc', row='optim_name', hue='Param_type', data=df_param_plot,
                    row_order=boxplot_order, kind='bar', size=3, aspect=4, margin_titles=True, legend=False, dodge=False) 
#    xcoords = [35.5,66.5]
    xcoords = [41.5,72.5]
    for xc in xcoords:
        g = g.map(plt.axvline, x=xc, ls='--', c='black') # add vertical lines to seperate parameter groups
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#    g.set_titles(row_template = '{row_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(hspace=0.5, wspace=0.2)
    plt.subplots_adjust(top=0.96,left=0.06,bottom=0.21,right=0.92)
    
    # entire figure title
#    g.fig.suptitle('Parameter uncertainty reduction for all optimizations with all parameters')   
    # one ylabel
    g.fig.text(x=0.005, y=0.6, verticalalignment='center', s='Parameter uncertainty reduction (fraction)', size=14, rotation=90)
    g.fig.text(x=0.25, y=0.97, horizontalalignment='center', s='Phenology', size=14, rotation=0)
    g.fig.text(x=0.625, y=0.97, horizontalalignment='center', s='Photosynthesis', size=14, rotation=0)
    g.fig.text(x=0.8, y=0.97, horizontalalignment='center', s='Post GPP', size=14, rotation=0)
    
    for nx, ax in enumerate(g.axes.flat):
        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=7) # set_xticklabels
        # This contains the right ylabel text
        txt = ax.texts[0]
        ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1],
                df_param_group.opt_param[nx+1],transform=ax.transAxes,va='center',fontsize='large')
    
    g.fig.savefig(outdir + 'Parameters/' + 'Parameters_uncertainty_reduction.png')
    plt.close()
    
    
    ## (NEE opt only) plotting all parameter uncertainty reductions
    ## ---------------------
    ## ---------------------
    df_param_plot = df_param_plot.loc[df_param_plot['optim'].isin([8,9,10,11,12,13,14])]
    df_param_plot = df_param_plot.drop('Param_type', 1)
    df_param_plot = pd.merge(df_param_plot, param_list, how='right', on='Param_name')
    df_param_plot = df_param_plot.sort_values(['Param_type','Param_name'])
    
    df_param_nee_uncer = df_param_plot.copy() # mark the parameters for NEE opt with all param
    
    #plot parameter deviation bar plot
    sns.set(font_scale=1, style='white')
    if any('paramP7' in s for s in optim_name): 
        g = sns.catplot(x='Param_name', y='uncer_reduc', row='optim', hue='optim_name', data=df_param_plot,
                    kind='bar', height=0.5, aspect=8, palette=['grey'], margin_titles=True, legend=False, dodge=False) 
    else:
        g = sns.catplot(x='Param_name', y='uncer_reduc', row='optim_name', hue='optim_name', data=df_param_plot,
                    kind='bar', height=0.5, aspect=8, palette=['grey'], margin_titles=True, legend=False, dodge=False) 
#    xcoords = [35.5,66.5]
    xcoords = [41.5,72.5]
    for xc in xcoords:
        g = g.map(plt.axvline, x=xc, ls='--', c='black') # add vertical lines to seperate parameter groups
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#    g.set_titles(row_template = '{row_name}')
    g.set(xlabel='', ylabel='')
    g.fig.set_figwidth(20)
    g.fig.set_figheight(5)
    plt.subplots_adjust(hspace=0.65, wspace=0.2)
    plt.subplots_adjust(top=0.98,left=0.05,bottom=0.62,right=0.96)
    
#    # entire figure title
#    g.fig.suptitle('Parameter uncertainty reduction for NEE optimizations with all parameter sets')   
#    # one ylabel
    g.fig.text(x=0.005, y=0.6, verticalalignment='center', s='Parameter uncertainty reduction', size=16, rotation=90)
#    g.fig.text(x=0.25, y=0.92, horizontalalignment='center', s='Phenology', size=14, rotation=0)
#    g.fig.text(x=0.64, y=0.92, horizontalalignment='center', s='Photosynthesis', size=14, rotation=0)
#    g.fig.text(x=0.89, y=0.92, horizontalalignment='center', s='Post GPP', size=14, rotation=0)
    
    for nx, ax in enumerate(g.axes.flat):
        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=12) # set_xticklabels
        # This contains the right ylabel text
        txt = ax.texts[0]
        ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1],
                df_param_group.param[nx+1],transform=ax.transAxes,va='center',fontsize='large')
    
    g.fig.savefig(outdir + 'Parameters/' + 'Parameters_uncertainty_reduction_NEE_opt.png')
    plt.close()
    
    
    
    #plotting all parameter deviations with asterisk to mark the uncertainty reduced parameters
    ## ---------------------
    ## ---------------------
    df_param_plot = df_param_all.loc[df_param_all['Paramtype'] == 'deviation']
    df_param_plot = df_param_plot.drop('Paramtype', 1)
    df_param_plot = df_param_plot.drop('Param_ID', 1)
    df_param_plot = df_param_plot.drop('Site', 1)
    
    df_param_plot['median_deviation'] = df_param_plot.groupby(['Param_name','optim']).Paramvalue.transform('median')
    df_param_plot = df_param_plot.drop_duplicates(subset='median_deviation')
        
    # read param list to add the param types
    param_list = pd.read_excel(param_file)
    df_optim = pd.DataFrame({'optim':np.repeat([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],len(param_list))})
    param_list = param_list.append([param_list]*19,ignore_index=True)
    param_list = pd.concat([param_list,df_optim], axis=1)
    
    # merge the data frames to have param type
    df_param_plot = pd.merge(df_param_plot, param_list, how='right', on=['Param_name','optim'])
    df_param_plot = df_param_plot.sort_values(['optim','Param_type','Param_name'], ascending=True)
    
    
    df_param_plot = df_param_plot.drop('Paramvalue', 1)
    
    df_param_yt = df_param_plot.copy() # mark the parameters for NEE opt with all param
    
    df_param_plot = pd.merge(df_param_plot, df_param_xt, how='left', on=['Param_name','optim','optim_name','Param_type'])
#    df_param_plot = pd.merge(df_param_plot, df_param_xt, left_on=['Param_name','optim','optim_name','Param_type'],
#                             right_on=['Param_name','optim','optim_name','Param_type'],how='outer')
#    df_param_plot.iloc[70:90,:]
    
    threshold = 50
    df_param_plot['asterisk'] = 1
    df_param_plot.loc[(df_param_plot['uncer_reduc'] > threshold), 'asterisk'] = 0.65
#    df_param_plot.loc[(df_param_plot['uncer_reduc'] > threshold), 
#           'asterisk'] = df_param_plot.loc[(df_param_plot['uncer_reduc'] > threshold), 
#           'median_deviation']*1.2
    df_param_plot['size'] = 0
    df_param_plot.loc[(df_param_plot['asterisk'] == 0.65),'size'] = 2
       
#    df_param_plot = df_param_plot.loc[df_param_plot['optim'].isin([1,2])]
#    df_param_plot = df_param_plot.reset_index()
#    df_param_plot.iloc[0:50,:]
#    df_param_plot.loc[:, df_param_plot.isna().any()]
#    #plot parameter deviation bar plot
#    g = sns.FacetGrid(df_param_plot, row='optim', hue='optim_name', height=12, aspect=2) 
#    g = g.map(sns.barplot, 'Param_name', 'median_deviation', data=df_param_plot)
#    g = g.map(sns.scatterplot, 'Param_name', 'asterisk', marker='*', size='size', data=df_param_plot)
    
    sns.set(font_scale=1, style='white')
    if any('paramP7' in s for s in optim_name): 
        g = sns.catplot(x='Param_name', y='median_deviation', row='optim', hue='optim_name', data=df_param_plot,
                    kind='bar', size=3, aspect=4, margin_titles=True, legend=False, dodge=False) 
    else: 
        g = sns.catplot(x='Param_name', y='median_deviation', row='optim_name', hue='Param_type', data=df_param_plot,
                    row_order=boxplot_order, kind='bar', size=3, aspect=4, margin_titles=True, legend=False, dodge=False) 
    xcoords = [41.5,72.5]
    for xc in xcoords:
        g = g.map(plt.axvline, x=xc, ls='--', c='black') # add vertical lines to seperate parameter groups
#    g = g.map(plt.scatter, x='Param_name', y='asterisk', data=df_param_plot, marker='*', s=2)
    g.map(sns.scatterplot, 'Param_name', 'asterisk', marker='*', size='size', data=df_param_plot)
    
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#    g.set_titles(row_template = '{row_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(hspace=0.4, wspace=0.2)
    plt.subplots_adjust(top=0.96,left=0.065,bottom=0.25,right=0.92)
    
#    # modify the legends
#    plt.legend(title='Optimization', bbox_to_anchor=(0.8, 2), ncol=2, fontsize=8)
#    new_labels = ['GPP+$R_{eco}$','NEE','GPP','$R_{eco}$']
#    for ax in g.axes.flat:
#        leg = ax.get_legend()
#        if not leg is None: break
#    # or legend may be on a figure
#    if leg is None: leg = ax.legend
#    
#    for t, l in zip(leg.texts, new_labels): t.set_text(l)
    
    # entire figure title
#    g.fig.suptitle('Parameter deviations for all optimizations with all parameters ( * indicates uncertainty reduction > 50%)')   
    # one ylabel
    g.fig.text(x=0.005, y=0.6, verticalalignment='center', s='Parameter deviations [(post-prior)/(max-min)]', size=14, rotation=90)
    g.fig.text(x=0.27, y=0.97, horizontalalignment='center', s='Phenology', size=14, rotation=0)
    g.fig.text(x=0.63, y=0.97, horizontalalignment='center', s='Photosynthesis', size=14, rotation=0)
    g.fig.text(x=0.85, y=0.97, horizontalalignment='center', s='Post GPP', size=14, rotation=0)
    g.set(ylim=(-0.5, 0.75), yticks=[-0.5,0.5])
    g.set(xlim=(-0.5, 86.5))
    
    for nx, ax in enumerate(g.axes.flat):
        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=8.5) # set_xticklabels
        
        # This contains the right ylabel text
        txt = ax.texts[0]
        ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1],
                df_param_group.opt_param[nx+1],transform=ax.transAxes,va='center',fontsize='large')
    
    
    g.fig.savefig(outdir + 'Parameters/' + 'Parameters_deviation_with_uncertainty_reduction.png')
    plt.close()
    
    
    # (NEE opt only) plotting all parameter deviations with asterisk to mark the uncertainty reduced parameters
    ## ---------------------
    ## ---------------------
    df_param_plot_nee = df_param_plot.loc[df_param_plot['optim'].isin([8,9,10,11,12,13,14])]
    
    sns.set(font_scale=1, style='white')
    if any('paramP7' in s for s in optim_name): 
        g = sns.catplot(x='Param_name', y='median_deviation', row='optim', hue='optim_name', data=df_param_plot_nee,
                    kind='bar', height=0.5, aspect=8, margin_titles=True, legend=False, dodge=False) 
    else: 
        g = sns.catplot(x='Param_name', y='median_deviation', row='optim_name', hue='Param_type', data=df_param_plot_nee,
                    kind='bar', height=0.5, aspect=8, margin_titles=True, legend=False, dodge=False) 
    xcoords = [41.5,72.5]
    for xc in xcoords:
        g = g.map(plt.axvline, x=xc, ls='--', c='black') # add vertical lines to seperate parameter groups
#    g = g.map(plt.scatter, x='Param_name', y='asterisk', data=df_param_plot_nee, marker='*', s=2)
    g = g.map(sns.scatterplot, 'Param_name', 'asterisk', marker='*', size='size', data=df_param_plot_nee)
    
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#    g.set_titles(row_template = '{row_name}')
    g.set(xlabel='', ylabel='')
    g.fig.set_figwidth(20)
    g.fig.set_figheight(5)
    plt.subplots_adjust(hspace=0.3, wspace=0.2)
    plt.subplots_adjust(top=0.92,left=0.05,bottom=0.02,right=0.96)
    
    # entire figure title
#    g.fig.suptitle('Parameter deviations for NEE optimizations with all parameter sets ( * indicates uncertainty reduction > 50%)')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='Parameter deviations [(post-prior)/(max-min)]', size=16, rotation=90)
    g.fig.text(x=0.25, y=0.94, horizontalalignment='center', s='Phenology', size=16, rotation=0)
    g.fig.text(x=0.65, y=0.94, horizontalalignment='center', s='Photosynthesis', size=16, rotation=0)
    g.fig.text(x=0.90, y=0.94, horizontalalignment='center', s='Post C uptake', size=16, rotation=0)
    g.set(ylim=(-0.5, 0.75), yticks=[-0.5,0.5])
    g.set(xlim=(-0.5, 86.5))
    g.set(xticklabels=[])
    
    for nx, ax in enumerate(g.axes.flat):
#        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=7) # set_xticklabels
        
        # This contains the right ylabel text
        txt = ax.texts[0]
        ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1],
                df_param_group.param[nx+1],transform=ax.transAxes,va='center',fontsize='large')
    
    
    g.fig.savefig(outdir + 'Parameters/' + 'Parameters_deviation_with_uncertainty_reduction_NEE_opt.png')
    plt.close()
    

    ## (NEE opt only) plotting all parameter uncertainty reductions
    #(but only the parameters with non-zero deviations)
    ## ---------------------
    ## ---------------------
    df_param_plot_nee_uncer = pd.merge(df_param_nee_uncer, df_param_yt, how='left', on=['Param_name','optim','optim_name','Param_type'])
    df_param_plot_nee_uncer = (df_param_plot_nee_uncer[~df_param_plot_nee_uncer.median_deviation.isnull()])
    
    #plot parameter deviation bar plot
    sns.set(font_scale=1, style='white')
    if any('paramP7' in s for s in optim_name): 
        g = sns.catplot(x='Param_name', y='uncer_reduc', row='optim', hue='optim_name', data=df_param_plot_nee_uncer,
                    kind='bar', height=0.5, aspect=8, palette=['grey'], margin_titles=True, legend=False, dodge=False) 
    else:
        g = sns.catplot(x='Param_name', y='uncer_reduc', row='optim_name', hue='optim_name', data=df_param_plot_nee_uncer,
                    kind='bar', height=0.5, aspect=8, palette=['grey'], margin_titles=True, legend=False, dodge=False) 
#    xcoords = [35.5,66.5]
    xcoords = [35.5,66.5]
    for xc in xcoords:
        g = g.map(plt.axvline, x=xc, ls='--', c='black') # add vertical lines to seperate parameter groups
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#    g.set_titles(row_template = '{row_name}')
    g.set(xlabel='', ylabel='')
    g.fig.set_figwidth(20)
    g.fig.set_figheight(5)
    plt.subplots_adjust(hspace=0.8, wspace=0.2)
    plt.subplots_adjust(top=0.98,left=0.05,bottom=0.58,right=0.96)
    
#    # entire figure title
#    g.fig.suptitle('Parameter uncertainty reduction for NEE optimizations with all parameter sets')   
#    # one ylabel
    g.fig.text(x=0.005, y=0.6, verticalalignment='center', s='Parameter uncertainty reduction', size=16, rotation=90)
#    g.fig.text(x=0.25, y=0.92, horizontalalignment='center', s='Phenology', size=14, rotation=0)
#    g.fig.text(x=0.64, y=0.92, horizontalalignment='center', s='Photosynthesis', size=14, rotation=0)
#    g.fig.text(x=0.89, y=0.92, horizontalalignment='center', s='Post GPP', size=14, rotation=0)
    
    for nx, ax in enumerate(g.axes.flat):
        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=11) # set_xticklabels
        # This contains the right ylabel text
        txt = ax.texts[0]
        ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1],
                df_param_group.param[nx+1],transform=ax.transAxes,va='center',fontsize='large')
    
    g.fig.savefig(outdir + 'Parameters/' + 'Parameters_uncertainty_reduction_NEE_opt_subset.png')
    plt.close()
    
    
    #plotting all parameter deviations with asterisk to mark the uncertainty reduced parameters 
    #(but only the parameters with non-zero deviations)
    ## ---------------------
    ## ---------------------
    df_param_plot = (df_param_plot[~df_param_plot.median_deviation.isnull()])
    sns.set(font_scale=1, style='white')
    if any('paramP7' in s for s in optim_name): 
        g = sns.catplot(x='Param_name', y='median_deviation', row='optim', hue='optim_name', data=df_param_plot,
                    kind='bar', size=3, aspect=4, margin_titles=True, legend=False, dodge=False) 
    else: 
        g = sns.catplot(x='Param_name', y='median_deviation', row='optim_name', hue='Param_type', data=df_param_plot,
                    row_order=boxplot_order, kind='bar', size=3, aspect=4, margin_titles=True, legend=False, dodge=False) 
    xcoords = [35.5,66.5]
    for xc in xcoords:
        g = g.map(plt.axvline, x=xc, ls='--', c='black') # add vertical lines to seperate parameter groups
#    g = g.map(plt.scatter, x='Param_name', y='asterisk', data=df_param_plot, marker='*', s=2)
    g.map(sns.scatterplot, 'Param_name', 'asterisk', marker='*', size='size', data=df_param_plot)
    
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#    g.set_titles(row_template = '{row_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(hspace=0.4, wspace=0.2)
    plt.subplots_adjust(top=0.96,left=0.065,bottom=0.25,right=0.92)
    
    # one ylabel
    g.fig.text(x=0.005, y=0.6, verticalalignment='center', s='Parameter deviations [(post-prior)/(max-min)]', size=14, rotation=90)
    g.fig.text(x=0.27, y=0.97, horizontalalignment='center', s='Phenology', size=14, rotation=0)
    g.fig.text(x=0.65, y=0.97, horizontalalignment='center', s='Photosynthesis', size=14, rotation=0)
    g.fig.text(x=0.875, y=0.97, horizontalalignment='center', s='Post GPP', size=14, rotation=0)
    g.set(ylim=(-0.5, 0.75), yticks=[-0.5,0.5])
    g.set(xlim=(-0.5, 80.5))
    
    for nx, ax in enumerate(g.axes.flat):
        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=8.5) # set_xticklabels
        
        # This contains the right ylabel text
        txt = ax.texts[0]
        ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1],
                df_param_group.opt_param[nx+1],transform=ax.transAxes,va='center',fontsize='large')
    
    
    g.fig.savefig(outdir + 'Parameters/' + 'Parameters_deviation_with_uncertainty_reduction_subset.png')
    plt.close()
    
    
    # (NEE opt only) plotting all parameter deviations with asterisk to mark the uncertainty reduced parameters
    #(but only the parameters with non-zero deviations)
    ## ---------------------
    ## ---------------------
    df_param_plot_nee = (df_param_plot_nee[~df_param_plot_nee.median_deviation.isnull()])
    sns.set(font_scale=1, style='white')
    if any('paramP7' in s for s in optim_name): 
        g = sns.catplot(x='Param_name', y='median_deviation', row='optim', hue='optim_name', data=df_param_plot_nee,
                    kind='bar', height=0.5, aspect=8, margin_titles=True, legend=False, dodge=False) 
    else: 
        g = sns.catplot(x='Param_name', y='median_deviation', row='optim_name', hue='Param_type', data=df_param_plot_nee,
                    kind='bar', height=0.5, aspect=8, margin_titles=True, legend=False, dodge=False) 
    xcoords = [35.5,66.5]
    for xc in xcoords:
        g = g.map(plt.axvline, x=xc, ls='--', c='black') # add vertical lines to seperate parameter groups
#    g = g.map(plt.scatter, x='Param_name', y='asterisk', data=df_param_plot_nee, marker='*', s=2)
    g = g.map(sns.scatterplot, 'Param_name', 'asterisk', marker='*', size='size', data=df_param_plot_nee)
    
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#    g.set_titles(row_template = '{row_name}')
    g.set(xlabel='', ylabel='')
    g.fig.set_figwidth(15)
    g.fig.set_figheight(4)
    plt.subplots_adjust(hspace=0.5, wspace=0.2)
    plt.subplots_adjust(top=0.92,left=0.05,bottom=0.02,right=0.96)
    
    # entire figure title
#    g.fig.suptitle('Parameter deviations for NEE optimizations with all parameter sets ( * indicates uncertainty reduction > 50%)')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='Parameter deviations [(post-prior)/(max-min)]', size=14, rotation=90)
    g.fig.text(x=0.25, y=0.94, horizontalalignment='center', s='Phenology', size=14, rotation=0)
    g.fig.text(x=0.65, y=0.94, horizontalalignment='center', s='Photosynthesis', size=14, rotation=0)
    g.fig.text(x=0.9, y=0.94, horizontalalignment='center', s='Post GPP', size=14, rotation=0)
    g.set(ylim=(-0.5, 0.75), yticks=[-0.5,0.5])
    g.set(xlim=(-0.5, 80.5))
    g.set(xticklabels=[])
    
    for nx, ax in enumerate(g.axes.flat):
#        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=7) # set_xticklabels
        
        # This contains the right ylabel text
        txt = ax.texts[0]
        ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1],
                df_param_group.param[nx+1],transform=ax.transAxes,va='center',fontsize='large')
    
    
    g.fig.savefig(outdir + 'Parameters/' + 'Parameters_deviation_with_uncertainty_reduction_NEE_opt_subset.png')
    plt.close()    
    

#%%    
    ## plotting the significant parameter deviations with a thresholding
    ## ---------------------
    ## ---------------------
    df_param_plot = df_param_all.loc[df_param_all['Paramtype'] == 'deviation']
    df_param_plot = df_param_plot.drop('Paramtype', 1)
    df_param_plot = df_param_plot.drop('Param_ID', 1)
    df_param_plot = df_param_plot.drop('Site', 1)
    
    df_param_plot['median_deviation'] = df_param_plot.groupby(['Param_name','optim']).Paramvalue.transform('median')
    df_param_plot = df_param_plot.drop_duplicates(subset='median_deviation')
    
    # read param list to add the param types
    param_list = pd.read_excel(param_file)
    
    # merge the data frames to have param type
    df_param_plot = pd.merge(df_param_plot, param_list, how='left', on='Param_name')
    df_param_plot = df_param_plot.sort_values(['Param_type','Param_name'], ascending=True)
    
    df_param_sub = df_param_plot.drop('Paramvalue', 1)
    
    threshold = 0.15
    df_param_sub['deviation_thresh'] = np.nan
    df_param_sub.loc[(df_param_sub['median_deviation'] > threshold) | 
           (df_param_sub['median_deviation'] < -threshold), 
           'deviation_thresh'] = df_param_sub.loc[(df_param_sub['median_deviation'] > threshold) | 
           (df_param_sub['median_deviation'] < -threshold), 
           'median_deviation']
            
    #df_param_sub[df_param_sub['deviation_thresh'].notna()]
    
    #plot parameter deviation bar plot
    sns.set(font_scale=1, style='white')
    if any('paramP7' in s for s in optim_name): 
        g = sns.catplot(x='Param_name', y='deviation_thresh', row='optim', data=df_param_sub,
                    kind='bar', height=12, aspect=2, margin_titles=True, legend=False, dodge=False) 
    else:
        g = sns.catplot(x='Param_name', y='deviation_thresh', row='optim_name', data=df_param_sub,
                    kind='bar', height=12, aspect=2, margin_titles=True, legend=False, dodge=False) 
    xcoords = [35.5,66.5]
    for xc in xcoords:
        g = g.map(plt.axvline, x=xc, ls='--', c='black') # add a horizontal line at ratio = 1.0
    
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#    g.set_titles(row_template = '{row_name}')
    g.set(xlabel='', ylabel='')
    # entire figure title
    plt.subplots_adjust(hspace=0.2, wspace=0.2)
    plt.subplots_adjust(top=0.94,left=0.05,bottom=0.22,right=0.96)
    
    g.fig.suptitle('Dominant parameters with a threshold: Parameter deviation < -' + str(threshold) + ' & > ' + str(threshold))
    # one ylabel
    g.fig.text(x=0.005, y=0.6, verticalalignment='center', s='Parameter deviations [(post-prior)/(max-min)]', size=12, rotation=90)
    g.fig.text(x=0.25, y=0.94, horizontalalignment='center', s='Phenology', size=12, rotation=0)
    g.fig.text(x=0.65, y=0.94, horizontalalignment='center', s='Photosynthesis', size=12, rotation=0)
    g.fig.text(x=0.90, y=0.94, horizontalalignment='center', s='Post GPP', size=12, rotation=0)
    
#    df_param_xt = df_param_sub.loc[df_param_sub['optim'] == 1] # mark the parameters for NEE opt with all param
#    df_param_xt = df_param_xt.reset_index()
#    df_param_xt = df_param_xt[['Param_name','deviation_thresh']].dropna()
    
    #show all xlabels with Dominant parameters marked as ------------------
    for nx, ax in enumerate(g.axes.flat):
        # This contains the right ylabel text
        txt = ax.texts[0]
        ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1],
                df_param_group.param[nx+1],transform=ax.transAxes,va='center',fontsize='large')
        
        labels = ax.get_xticklabels() # get x labels
#        for i,l in enumerate(labels):
#            if(i in df_param_xt.index): labels[i] = '--------------- '+labels[i].get_text() # skip even labels
        ax.set_xticklabels(labels, rotation=90, fontsize=7) # set new labels 
    
    g.fig.savefig(outdir + 'Parameters/' + 'Parameters_deviation_subset.png')
    plt.close()

    
    
    ## plotting the site-wise parameters with range, error bars
    # choose the parameters for PFT=11 (C4G=85%) for the US-Wjs site
    if 'US-Wjs' in df_param_all.Site.unique():
            param_wjs = df_param_all[df_param_all.Site == 'US-Wjs']
            param_wjs = param_wjs[~param_wjs.Param_ID.str.contains('__04')] # remove PFT=04 (15%) for US-Wjs
            df_param_all = df_param_all[~df_param_all.Site.str.contains('Wjs')]
            df_param_all = pd.concat([df_param_all,param_wjs], ignore_index=True)
        
        # choose the parameters for PFT=11 (C4G=44%) for the US-SRG site
    if 'US-SRG' in df_param_all.Site.unique():
            param_SRG = df_param_all[df_param_all.Site == 'US-SRG']
            param_SRG = param_SRG[~param_SRG.Param_ID.str.contains('__06')] # remove PFT=06 (11%) for US-SRG
            df_param_all = df_param_all[~df_param_all.Site.str.contains('SRG')]
            df_param_all = pd.concat([df_param_all,param_SRG], ignore_index=True)    
        
        # choose the parameters for PFT=06 (TeBD=35%) for the US-SRM site
    if 'US-SRM' in df_param_all.Site.unique():
            param_SRM = df_param_all[df_param_all.Site == 'US-SRM']
            param_SRM = param_SRM[~param_SRM.Param_ID.str.contains('__11')] # remove PFT=11 (15%) for US-SRM
            df_param_all = df_param_all[~df_param_all.Site.str.contains('SRM')]
            df_param_all = pd.concat([df_param_all,param_SRM], ignore_index=True)
        
        # choose the parameters for PFT=05 (TeBE=55%) for the US-Ses site
    if 'US-Ses' in df_param_all.Site.unique():
            param_Ses = df_param_all[df_param_all.Site == 'US-Ses']
            param_Ses = param_Ses[~param_Ses.Param_ID.str.contains('__11')] # remove PFT=11 (25%) for US-Ses
            df_param_all = df_param_all[~df_param_all.Site.str.contains('Ses')]
            df_param_all = pd.concat([df_param_all,param_Ses], ignore_index=True)
        
        # choose the parameters for PFT=04 (TeNE=60%) for the US-SRM site
    if 'US-Mpj' in df_param_all.Site.unique():
            param_Mpj = df_param_all[df_param_all.Site == 'US-Mpj']
            param_Mpj = param_Mpj[~param_Mpj.Param_ID.str.contains('__11')] # remove PFT=11 (20%) for US-Mpj
            df_param_all = df_param_all[~df_param_all.Site.str.contains('Mpj')]
            df_param_all = pd.concat([df_param_all,param_Mpj], ignore_index=True)
        
        #subset the dataframe
    df_param_plot = df_param_all[df_param_all['Paramtype'].isin(['Param_post','Param_default','Param_min','Param_max'])]
        
        # read param list to add the param types
    param_list = pd.read_excel(param_file)
    df_optim = pd.DataFrame({'Site':np.repeat(site_names,len(param_list))})
    param_list = param_list.append([param_list]*(len(site_names)-1),ignore_index=True)
    param_list = pd.concat([param_list,df_optim], axis=1)
        
        # merge the data frames to have param type
    df_param_plot = pd.merge(df_param_plot, param_list, how='right', on=['Param_name','Site'])
    df_param_plot = df_param_plot.sort_values(['Param_type','Param_name'], ascending=True)
    
#    df_param_group['opt_param'] = df_param_group["opt"] + df_param_group["param"]
    df_param_group.opt_param[0] = 'Prior'
    
    #plot all parameters
    for ns, ss in enumerate(site_names):
#        ns = 0
        ss = site_names[ns]
        df_param_plot_sub = df_param_plot.loc[df_param_plot['Site'] == ss]
        
        sns.set(font_scale=0.8, style='white')
        g = sns.relplot(x='optim', y='Paramvalue', hue='Paramtype', style='Paramtype', 
                size = 'Paramtype', col='Param_name',
                data=df_param_plot_sub, palette=['purple','grey','orange','orange'], 
                col_wrap=11, height=10, aspect=1,
                markers=["o","X",'v','^'], sizes=(10,30),
                facet_kws=dict(sharey=False),
                size_order = ['Param_post','Param_default','Param_min','Param_max'])  

        plt.xlim(-0.5, 21)
        plt.subplots_adjust(hspace=0.3, wspace=0.45)
        [plt.setp(ax.texts, text="") for ax in g.axes.flat]
        g.set_titles(col_template = '{col_name}')
        g.set(xlabel='', ylabel='')
        plt.subplots_adjust(top=0.97,left=0.035,bottom=0.07,right=0.99)
        g._legend.set_bbox_to_anchor([0.99,0.08])
        # replace labels
        new_labels = ['','Posterior parameters','Prior parameters','Parameter bounds','Parameter bounds']
        for t, l in zip(g._legend.texts, new_labels): t.set_text(l)
                    
#        g.fig.suptitle('Parameters ' + r"$\bf{FOR}$ " + ' site ' + site_names[ns] + r" $\bf{WITH}$ " + 'all optimizations')   
        g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='Parameter values', size=12, rotation=90)
        
        g.set(xticks=list(range(1, 21)))
        g.set_xticklabels(df_param_group.opt_param[1:], rotation=90, fontsize=4)
        
        
        g.fig.savefig(outdir + 'Parameters/' + 'Parameters_all_optim_' + site_names[ns] + '_Csink.png', dpi = 1000)
        plt.close()
        
    df_param_group  = df_param_group.drop('opt', 1)
    df_param_group  = df_param_group.drop('opt_param', 1)
        
#%% plotting Box-whisker plots for correlation coefficient and slope of regression comparison for all optimizations with all parameter groups
## ---------------------
## ---------------------
if corr_plot == True: 
    boxplot_order=[0,8,9,10,11,12,13,14,1,2,3,4,5,6,7,15,16,17,18,19,20]
    boxplot_order_sub=[0,8,9,10,11,12,13,14]
    pal2=['blue','orange','green']
    #assign site type (source / pivot / sink)
    flux_msd_all['site_type'] = 'Sink'
#    flux_msd_all.loc[flux_msd_all['Site'] == 'US-SO2', 'site_type'] = 'Pivot'
    flux_msd_all.loc[flux_msd_all['Site'] == 'US-Wkg', 'site_type'] = 'Pivot'
    flux_msd_all.loc[flux_msd_all['Site'] == 'US-Whs', 'site_type'] = 'Pivot'
    flux_msd_all.loc[flux_msd_all['Site'] == 'US-SRG', 'site_type'] = 'Pivot'
    flux_msd_all.loc[flux_msd_all['Site'] == 'US-SRM', 'site_type'] = 'Pivot'
    flux_msd_all.loc[flux_msd_all['Site'] == 'US-Seg', 'site_type'] = 'Pivot'
    
#    flux_msd_all.loc[flux_msd_all['Site'] == 'US-Seg', 'site_type'] = 'Source'
    flux_msd_all.loc[flux_msd_all['Site'] == 'US-Aud', 'site_type'] = 'Source'
    sns.set(font_scale=1.1, style='white')
    
    for nf, sf in enumerate(flux_names):
#        nf = 0
        sf = flux_names[nf]
        
        #Box-whisker plots for correlation coefficient
        flux_corr_all_plot = flux_corr_all.loc[flux_corr_all['Flux'] == sf]
        flux_corr_all_plot.optim[flux_corr_all_plot.Simulation == 'Prior'] = 0
        flux_corr_all_plot.optim_name[flux_corr_all_plot.Simulation == 'Prior'] = 'Prior'
        
        flux_corr_all_plot = pd.merge(flux_corr_all_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
        
        g = sns.boxplot(x='optim', y='Anomaly', hue='optim_name', data=flux_corr_all_plot, 
                        order=boxplot_order, hue_order=['Prior','nee','gpp_reco','gpp','reco'], dodge=False)  
        g.figure.set_size_inches(5,4)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.93,bottom=0.12,left=0.15,right=0.99)
        plt.legend(title='Optimization', bbox_to_anchor=(0.45, 0.25), ncol=2, fontsize=8)
        plt.suptitle(subplot_title[nf]+' '+sf)   
        plt.text(x=-3.9, y=(max(flux_corr_all_plot.Anomaly)+min(flux_corr_all_plot.Anomaly))/2, verticalalignment='center', s='Correlation coefficient (R)', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param, rotation=90, fontsize=10)
        new_labels = ['Prior','NEE','GPP+$R_{eco}$','GPP','$R_{eco}$']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)
        
        if nf > 0: 
            g.legend_.remove()

#        if sf == 'NEE':
#            plt.text(x=18, y=1.125, s='(a)', size=16) 
#        elif sf == 'GPP':
#            plt.text(x=18, y=1.09, s='(b)', size=16) 
#        else:
#            plt.text(x=18, y=1.09, s='(c)', size=16) 
        
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_anomaly_corrcoef_box_whisker.png')
        plt.close()
        
        #Box-whisker plots for correlation coefficient (only NEE opt)
        flux_corr_all_plot = flux_corr_all_plot[flux_corr_all_plot['optim_name'].isin(['Prior','nee'])]
        
        g = sns.boxplot(x='optim', y='Anomaly', hue='optim_name', data=flux_corr_all_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(5,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.93,bottom=0.12,left=0.17,right=0.99)
#        plt.subplots_adjust(top=0.93,bottom=0.12,left=0.27,right=0.99)
        plt.legend(title='Optimization', bbox_to_anchor=(0.8, 0.25), ncol=2, fontsize=14)
        plt.text(x=6.8, y=min(flux_corr_all_plot.Anomaly), s=subplot_title[nf], size=20)  
        plt.suptitle(sf)   
        if sf == 'NEE':
            plt.text(x=-2, y=(max(flux_corr_all_plot.Anomaly)+min(flux_corr_all_plot.Anomaly))/2, verticalalignment='center', s='Correlation coefficient (R)', size=16, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
        new_labels = ['Prior','Posterior']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)
        
        if nf > 0: 
            g.legend_.remove()

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_anomaly_corrcoef_sub_box_whisker.png')
        plt.close()

        
        ## ---------------------
        #Box-whisker plots for R^2
        flux_corr_all_plot = flux_corr_all.loc[flux_corr_all['Flux'] == sf]
        flux_corr_all_plot.optim[flux_corr_all_plot.Simulation == 'Prior'] = 0
        flux_corr_all_plot.optim_name[flux_corr_all_plot.Simulation == 'Prior'] = 'Prior'
        
        flux_corr_all_plot = pd.merge(flux_corr_all_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
        flux_corr_all_plot.Anomaly = flux_corr_all_plot.Anomaly ** 2
        
        g = sns.boxplot(x='optim', y='Anomaly', hue='optim_name', data=flux_corr_all_plot, 
                        order=boxplot_order, hue_order=['Prior','nee','gpp_reco','gpp','reco'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best', ncol=5)
        plt.suptitle(sf+' Flux')   
        plt.text(x=-1.6, y=max(flux_corr_all_plot.Anomaly)/2, verticalalignment='center', s='Correlation coefficient ($R^{2}$)', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param, rotation=90, fontsize=10)
        new_labels = ['Prior','NEE','GPP+$R_{eco}$','GPP','$R_{eco}$']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_anomaly_corrcoef^2_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for R^2 (only NEE opt)
        flux_corr_all_plot = flux_corr_all_plot[flux_corr_all_plot['optim_name'].isin(['Prior','nee'])]
        
        sns.set(font_scale=1.1, style='white')
        g = sns.boxplot(x='optim', y='Anomaly', hue='optim_name', data=flux_corr_all_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(5,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.99,bottom=0.12,left=0.13,right=0.99)
#        plt.legend(title='Optimization', bbox_to_anchor=(0.2,0.9), ncol=2)
#        plt.suptitle(sf+' Flux')   
        g.legend_.remove()
        plt.text(x=3.3, y=0.93, s='(e)', size=16)  
        plt.text(x=-1.65, y=max(flux_corr_all_plot.Anomaly)/2, verticalalignment='center', s='Correlation coefficient ($R^{2}$)', size=14, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
#        new_labels = ['Prior','NEE']
#        leg = g.axes.get_legend()
#        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_anomaly_corrcoef^2_sub_box_whisker.png')
        plt.close()
    
    
        ## ---------------------
        #Box-whisker plots for slope of regression of annual data
        flux_slope_all_plot = flux_slope_all.loc[flux_slope_all['Flux'] == sf]
        flux_slope_all_plot.optim[flux_slope_all_plot.Simulation == 'Prior'] = 0
        flux_slope_all_plot.optim_name[flux_slope_all_plot.Simulation == 'Prior'] = 'Prior'
        
        flux_slope_all_plot = pd.merge(flux_slope_all_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
        
        g = sns.boxplot(x='optim', y='Anomaly', hue='optim_name', data=flux_slope_all_plot, 
                        order=boxplot_order, hue_order=['Prior','nee','gpp_reco','gpp','reco'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(subplot_title[nf+3]+' '+sf)   
        plt.text(x=-1.6, y=(max(flux_slope_all_plot.Anomaly)+min(flux_slope_all_plot.Anomaly))/2, verticalalignment='center', s='Slope of linear least-squares regression', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param, rotation=90, fontsize=10)
        new_labels = ['Prior','NEE','GPP+$R_{eco}$','GPP','$R_{eco}$']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_anomaly_slope_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for slope of regression (only NEE opt)
        flux_slope_all_plot = flux_slope_all_plot[flux_slope_all_plot['optim_name'].isin(['Prior','nee'])]
        
        sns.set(font_scale=1.1, style='white')
        g = sns.boxplot(x='optim', y='Anomaly', hue='optim_name', data=flux_slope_all_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(5,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.93,bottom=0.12,left=0.17,right=0.99)
        g.legend_.remove()
#        plt.legend(title='Optimization', loc='best', fontsize=12)
#        plt.suptitle(subplot_title[nf+3])   
        plt.text(x=6.8, y=min(flux_slope_all_plot.Anomaly), s=subplot_title[nf+3], size=20)  
        if sf == 'NEE':
            plt.text(x=-2, y=(max(flux_slope_all_plot.Anomaly)+min(flux_slope_all_plot.Anomaly))/2, verticalalignment='center', s='Slope of linear least-square regression', size=16, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
#        new_labels = ['Prior','NEE']
#        leg = g.axes.get_legend()
#        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_anomaly_slope_sub_box_whisker.png')
        plt.close()
        
        
        ## ---------------------
        #Box-whisker plots for bias
        flux_msd_all_plot = flux_msd_all.loc[flux_msd_all['Flux'] == sf]
        
        flux_msd_all_plot.optim[flux_msd_all_plot.Sim == 'prior'] = 0
        flux_msd_all_plot.optim_name[flux_msd_all_plot.Sim == 'prior'] = 'Prior'
        
        flux_msd_all_plot = flux_msd_all_plot.drop('mse', 1)
        flux_msd_all_plot = flux_msd_all_plot.drop('r', 1)
        flux_msd_all_plot = pd.merge(flux_msd_all_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
        
        g = sns.boxplot(x='optim', y='bias', hue='optim_name', data=flux_msd_all_plot, 
                        order=boxplot_order, hue_order=['Prior','nee','gpp_reco','gpp','reco'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
#        plt.text(x=-1.6, y=max(flux_msd_all_plot.bias)/2, verticalalignment='center', s='$Bias^{2}$ contribution to MSD', size=12, rotation=90)
        plt.text(x=-1.6, y=max(flux_msd_all_plot.bias)/2, verticalalignment='center', s='$Bias^{2}$', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param, rotation=90, fontsize=10)
        new_labels = ['Prior','NEE','GPP+$R_{eco}$','GPP','$R_{eco}$']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_bias_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for variance
        g = sns.boxplot(x='optim', y='variance', hue='optim_name', data=flux_msd_all_plot, 
                        order=boxplot_order, hue_order=['Prior','nee','gpp_reco','gpp','reco'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
#        plt.text(x=-1.6, y=max(flux_msd_all_plot.variance)/2, verticalalignment='center', s='$Variance^{2}$ contribution to MSD', size=12, rotation=90)
        plt.text(x=-1.6, y=max(flux_msd_all_plot.variance)/2, verticalalignment='center', s='$Variance^{2}$', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param, rotation=90, fontsize=10)
        new_labels = ['Prior','NEE','GPP+$R_{eco}$','GPP','$R_{eco}$']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_variance_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for phase
        g = sns.boxplot(x='optim', y='phase', hue='optim_name', data=flux_msd_all_plot, 
                        order=boxplot_order, hue_order=['Prior','nee','gpp_reco','gpp','reco'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
        plt.text(x=-1.6, y=max(flux_msd_all_plot.phase)/2, verticalalignment='center', s='Phase contribution to MSD', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param, rotation=90, fontsize=10)
        new_labels = ['Prior','NEE','GPP+$R_{eco}$','GPP','$R_{eco}$']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_phase_box_whisker.png')
        plt.close()
        
        
        ## ---------------------
        #Box-whisker plots for bias, variance and phase side by side
        msd_melt = pd.melt(flux_msd_all_plot, id_vars=('Sim','Site','site_type','Flux','optim','optim_name','param'), var_name="variable", value_name="value")
        
        g = sns.boxplot(x='optim', y='value', hue='variable', data=msd_melt,order=boxplot_order,
                      hue_order=['bias','variance','phase'])  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        xcoords = [0.5,7.5,14.5,17.5]
        for xc in xcoords:
            plt.axvline(x=xc, ls='--', c='black') # add a horizontal line at ratio = 1.0
        
        g.text(x=0.08, y=-max(msd_melt.value)/20, horizontalalignment='center', s='Prior', size=11, rotation=0, weight='bold')
        g.text(x=3.5, y=-max(msd_melt.value)/20, horizontalalignment='center', s='NEE Optimization', size=11, rotation=0, weight='bold')
        g.text(x=11.5, y=-max(msd_melt.value)/20, horizontalalignment='center', s='GPP+$R_{eco}$ Optimization', size=11, rotation=0, weight='bold')
        g.text(x=16, y=-max(msd_melt.value)/20, horizontalalignment='center', s='GPP Optimization', size=11, rotation=0, weight='bold')
        g.text(x=19, y=-max(msd_melt.value)/20, horizontalalignment='center', s='$R_{eco}$ Optimization', size=11, rotation=0, weight='bold')
        
        plt.subplots_adjust(top=0.94,bottom=0.1,left=0.06,right=0.99)
        plt.legend(title='', loc='best')
        plt.suptitle(sf+' Flux')   
#        plt.text(x=-1.8, y=max(msd_melt.value)/2, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase contribution to MSD', size=12, rotation=90)
        plt.text(x=-1.8, y=max(msd_melt.value)/2, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase', size=12, rotation=90)
        plt.ylim(-max(msd_melt.value)/15,)
        
        g.set_xticklabels(df_param_group.param, rotation=90, fontsize=11)
        
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_bias_variance_phase_box_whisker.png')
        plt.close()
        
        
        ## ---------------------
        #Box-whisker plots for bias (only NEE opt)
        flux_msd_all_plot = flux_msd_all_plot[flux_msd_all_plot['optim_name'].isin(['Prior','nee'])]
        
        g = sns.boxplot(x='optim', y='bias', hue='optim_name', data=flux_msd_all_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
        plt.text(x=-0.95, y=max(flux_msd_all_plot.bias)/2, verticalalignment='center', s='$Bias^{2}$', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_bias_sub_box_whisker.png')
        plt.close()
        
        #Box-whisker plots for variance (only NEE opt)
        g = sns.boxplot(x='optim', y='variance', hue='optim_name', data=flux_msd_all_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
        plt.text(x=-0.95, y=max(flux_msd_all_plot.variance)/2, verticalalignment='center', s='$Variance^{2}$', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_variance_sub_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for phase (only NEE opt)
        g = sns.boxplot(x='optim', y='phase', hue='optim_name', data=flux_msd_all_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', bbox_to_anchor=(0.6,1))
        plt.suptitle(sf+' Flux')   
        plt.text(x=-0.95, y=max(flux_msd_all_plot.phase)/2, verticalalignment='center', s='Phase', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_phase_sub_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for bias, variance and phase side by side (only NEE opt)
        msd_melt = pd.melt(flux_msd_all_plot, id_vars=('Sim','Site','site_type','Flux','optim','optim_name','param'), var_name="variable", value_name="value")
        msd_melt = msd_melt.drop_duplicates(subset=['Sim','Site','site_type','Flux','optim','optim_name','param','variable'], keep="first")

        sns.set(font_scale=1.1, style='white')
        fig, axes = plt.subplots(3, 1, figsize=(4,15))
        axes = axes.flatten()
        g1 = sns.boxplot(x='optim', y='value', hue='variable', data=msd_melt[msd_melt['site_type'].isin(['Sink'])], dodge=True, ax=axes[0])  
        g1.set(xlabel='', ylabel='')
        g1.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g1.set(xticklabels=[])
        g1.legend_.remove()
#        g1.set(title = subplot_title[nf]+' '+sf)
#        g1.set(ylim=(-0.02*max(msd_melt.value), 0.02*max(msd_melt.value)+max(msd_melt.value)))
        g1.set(title = sf)
        if sf == 'NEE':
            g1.set(ylim=(-0.02*max(msd_melt.value), 0.02*max(msd_melt.value)+max(msd_melt.value)))
            g1.text(x=6.4, y=0.85, verticalalignment='center', s='(a)', size=12)        
        else: 
            g1.set(ylim=(-0.1,3.3))
        if sf == 'GPP': g1.text(x=6.6, y=3.1, verticalalignment='center', s='(b)', size=12)        
        if sf == 'Reco': g1.text(x=6.6, y=3.1, verticalalignment='center', s='(c)', size=12)        
        
        g2 = sns.boxplot(x='optim', y='value', hue='variable', data=msd_melt[msd_melt['site_type'].isin(['Pivot'])], dodge=True, ax=axes[1])  
        g2.set(xlabel='', ylabel='')
        g2.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g2.set(xticklabels=[])
        g2.legend_.remove()
#        g2.set(ylim=(-0.02*max(msd_melt.value), 0.02*max(msd_melt.value)+max(msd_melt.value)))
        if sf == 'NEE':
            g2.set(ylim=(-0.02*max(msd_melt.value), 0.02*max(msd_melt.value)+max(msd_melt.value)))
            g2.text(x=6.4, y=0.85, verticalalignment='center', s='(d)', size=12)        
        else: g2.set(ylim=(-0.1,3.3))
        if sf == 'GPP': g2.text(x=6.6, y=3.1, verticalalignment='center', s='(e)', size=12)        
        if sf == 'Reco': g2.text(x=6.6, y=3.1, verticalalignment='center', s='(f)', size=12)        
        
        g3 = sns.barplot(x="optim", y="value", hue="variable", data=msd_melt[msd_melt['site_type'].isin(['Source'])], dodge=True, ax=axes[2])
#        g3 = sns.stripplot(x="optim", y="value", edgecolor="none", hue="variable", linewidth=2,
#                  marker='_', data=msd_melt[msd_melt['site_type'].isin(['Source'])], dodge=True, ax=axes[2])
        g3.set(xlabel='', ylabel='')
        g3.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g3.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        g3.legend_.remove()
#        g3.set(ylim=(-0.02*max(msd_melt.value), 0.02*max(msd_melt.value)+max(msd_melt.value)))
        if sf == 'NEE':
            g3.set(ylim=(-0.02*max(msd_melt.value), 0.02*max(msd_melt.value)+max(msd_melt.value)))
            g3.text(x=6.4, y=0.85, verticalalignment='center', s='(g)', size=12) 
        else: 
            g3.set(ylim=(-0.1,3.3))
            
        if sf == 'GPP': g3.text(x=6.6, y=3.1, verticalalignment='center', s='(h)', size=12)        
        if sf == 'Reco': g3.text(x=6.6, y=3.1, verticalalignment='center', s='(i)', size=12)        
        
        if sf == 'NEE': 
            plt.subplots_adjust(top=0.95,bottom=0.07,left=0.2,right=0.99)
        elif sf == 'GPP': 
            plt.subplots_adjust(top=0.95,bottom=0.07,left=0.12,right=0.99)
        elif sf == 'Reco': 
            plt.subplots_adjust(top=0.95,bottom=0.07,left=0.12,right=0.93)
            
        plt.subplots_adjust(wspace=0, hspace=0.03)
        if sf == 'Reco':
            plt.legend(title='', loc='upper center', fontsize=12)
#            plt.subplots_adjust(top=0.95,bottom=0.07,left=0.15,right=0.93)
            plt.text(x=7.6, y=1.5, verticalalignment='center', s='Source', size=14, rotation=-90)        
            plt.text(x=7.6, y=5, verticalalignment='center', s='Pivot', size=14, rotation=-90)        
            plt.text(x=7.6, y=8.5, verticalalignment='center', s='Sink', size=14, rotation=-90)        
        if sf == 'NEE':
            plt.text(x=-2.3, y=1.4, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase', size=14, rotation=90)        
        
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_bias_variance_phase_sub_box_whisker.png')
        plt.close()
        
        if sf=='NEE':
            sns.set(font_scale=1.2, style='white')
            fig, axes = plt.subplots(1, 3, figsize=(15,5))
            axes = axes.flatten()
            g1 = sns.boxplot(x='optim', y='value', hue='variable', data=msd_melt[msd_melt['site_type'].isin(['Sink'])], dodge=True, ax=axes[0])  
            g1.set(xlabel='', ylabel='')
            g1.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
            g1.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
            g1.legend_.remove()
            g1.set(title = 'Sink')
            g1.set(ylim=(-0.02*max(msd_melt.value), 0.02*max(msd_melt.value)+max(msd_melt.value)))
            g1.text(x=6.4, y=0.85, verticalalignment='center', s='(a)', size=16)        
            g1.add_patch(Rectangle((-0.5, -0.015*max(msd_melt.value)), 2, 0.03*max(msd_melt.value)+max(msd_melt.value),linewidth=3,linestyle='dashed',edgecolor='grey',facecolor='none'))
            
            g2 = sns.boxplot(x='optim', y='value', hue='variable', data=msd_melt[msd_melt['site_type'].isin(['Pivot'])], dodge=True, ax=axes[1])  
            g2.set(xlabel='', ylabel='')
            g2.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
            g2.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
            g2.legend_.remove()
            g2.set(title = 'Pivot')
            g2.set(ylim=(-0.02*max(msd_melt.value), 0.02*max(msd_melt.value)+max(msd_melt.value)))
            g2.text(x=6.4, y=0.85, verticalalignment='center', s='(b)', size=16)        
            g2.add_patch(Rectangle((-0.5, -0.015*max(msd_melt.value)), 2, 0.03*max(msd_melt.value)+max(msd_melt.value),linewidth=3,linestyle='dashed',edgecolor='grey',facecolor='none'))
            
            g3 = sns.barplot(x="optim", y="value", hue="variable", data=msd_melt[msd_melt['site_type'].isin(['Source'])], dodge=True, ax=axes[2])
            g3.set(xlabel='', ylabel='')
            g3.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
            g3.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
            g3.set(title = 'Source')
            g3.set(ylim=(-0.02*max(msd_melt.value), 0.02*max(msd_melt.value)+max(msd_melt.value)))
            g3.text(x=6.4, y=0.85, verticalalignment='center', s='(c)', size=16)        
            g3.add_patch(Rectangle((-0.5, -0.015*max(msd_melt.value)), 2, 0.03*max(msd_melt.value)+max(msd_melt.value),linewidth=3,linestyle='dashed',edgecolor='grey',facecolor='none'))
            
            plt.subplots_adjust(top=0.94,bottom=0.13,left=0.05,right=0.99)
            plt.subplots_adjust(wspace=0.12, hspace=0)
            plt.legend(title='', loc='upper center', fontsize=12)
            plt.text(x=-19.7, y=0.4, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase', size=14, rotation=90)        
            
            plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_only_daily_bias_variance_phase_sub_box_whisker.png')
            plt.close()
        
        
        
#        Previous plot with 2 source sites and box plots
#        sns.set(font_scale=1.25, style='white')
#        g = sns.catplot(x='optim', y='value', hue="variable", row='site_type', 
#                data=msd_melt, kind="box", size=3, height=1, aspect=1.5, margin_titles=True, legend=False)
#        g.set(xlabel='', ylabel='')
#        g = g.map(plt.axvline, x=0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
#        [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#                
#        plt.subplots_adjust(top=0.95,bottom=0.07,left=0.16,right=0.99)
#        if sf == 'Reco':
#            plt.legend(title='', loc='upper center', fontsize=11)
#            g.set_titles(row_template = '{row_name}')
#            plt.subplots_adjust(top=0.95,bottom=0.07,left=0.16,right=0.93)
#        
#        plt.suptitle(subplot_title[nf]+' '+sf, fontsize=14)
#        
#        if sf == 'NEE':
#            plt.text(x=-1.9, y=3.5*max(msd_melt.value)/2, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase', size=14, rotation=90)        
#        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
#        
#        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_bias_variance_phase_sub_box_whisker.png')
#        plt.close()
        
        
        #Box-whisker plots for bias, variance and phase side by side (only NEE opt and Prior + P1)
        msd_melt_p1 = msd_melt[msd_melt['param'].isin(['Prior','P1'])]
        
#        msd_melt_p1 = msd_melt_p1.drop_duplicates(subset=['Sim','Site','site_type','Flux','optim','optim_name','param','variable'], keep="first")
#        msd_melt_p1_sub = msd_melt_p1[msd_melt_p1['Site'].isin(['US-Aud'])]
#        pd.set_option('display.max_rows', 500)
#        pd.set_option('display.max_columns', 500)
        
        sns.set(font_scale=1, style='white')
        fig, axes = plt.subplots(3, 1, figsize=(2,5))
        axes = axes.flatten()
        g1 = sns.boxplot(x='optim', y='value', hue='variable', data=msd_melt_p1[msd_melt_p1['site_type'].isin(['Sink'])], dodge=True, ax=axes[0])  
        g1.set(xlabel='', ylabel='')
        g1.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g1.set(xticklabels=[])
        g1.legend_.remove()
#        g1.set(title = subplot_title[nf]+' '+sf)
        g1.set(title = sf)
        if sf == 'NEE':
            g1.set(ylim=(-0.02*max(msd_melt_p1.value), 0.02*max(msd_melt_p1.value)+max(msd_melt_p1.value)))
            g1.text(x=1.15, y=0.75, verticalalignment='center', s='(a)', size=12)        
        else: 
            g1.set(ylim=(-0.06,2.88))
        if sf == 'GPP': g1.text(x=1.2, y=2.6, verticalalignment='center', s='(b)', size=12)        
        if sf == 'Reco': g1.text(x=1.2, y=2.6, verticalalignment='center', s='(c)', size=12)        
                
        g2 = sns.boxplot(x='optim', y='value', hue='variable', data=msd_melt_p1[msd_melt_p1['site_type'].isin(['Pivot'])], dodge=True, ax=axes[1])  
        g2.set(xlabel='', ylabel='')
        g2.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g2.set(xticklabels=[])
        g2.legend_.remove()
        if sf == 'NEE':
            g2.set(ylim=(-0.02*max(msd_melt_p1.value), 0.02*max(msd_melt_p1.value)+max(msd_melt_p1.value)))
            g2.text(x=1.15, y=0.75, verticalalignment='center', s='(d)', size=12)        
        else: g2.set(ylim=(-0.06,2.88))
        if sf == 'GPP': g2.text(x=1.2, y=2.6, verticalalignment='center', s='(e)', size=12)        
        if sf == 'Reco': g2.text(x=1.2, y=2.6, verticalalignment='center', s='(f)', size=12)        
        
        g3 = sns.barplot(x="optim", y="value", hue="variable", data=msd_melt_p1[msd_melt_p1['site_type'].isin(['Source'])], dodge=True, ax=axes[2])
#        g3 = sns.stripplot(x="optim", y="value", edgecolor="none", hue="variable", linewidth=2,
#                  marker='_', data=msd_melt_p1[msd_melt_p1['site_type'].isin(['Source'])], dodge=True, ax=axes[2])
        g3.set(xlabel='', ylabel='')
        g3.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g3.set_xticklabels(df_param_group.param[0:2], rotation=90, fontsize=12)
        g3.legend_.remove()
        if sf == 'NEE':
            g3.set(ylim=(-0.02*max(msd_melt_p1.value), 0.02*max(msd_melt_p1.value)+max(msd_melt_p1.value)))
            g3.text(x=1.15, y=0.75, verticalalignment='center', s='(g)', size=12)        
        else: g3.set(ylim=(-0.06,2.88))
        if sf == 'GPP': g3.text(x=1.2, y=2.6, verticalalignment='center', s='(h)', size=12)        
        if sf == 'Reco': g3.text(x=1.2, y=2.6, verticalalignment='center', s='(i)', size=12)        
        
        plt.subplots_adjust(wspace=0, hspace=0.03)
        if sf == 'GPP': plt.subplots_adjust(top=0.95,bottom=0.11,left=0.15,right=0.99)
        if sf == 'Reco':
            plt.legend(title='', loc='best', fontsize=10)
            plt.subplots_adjust(top=0.95,bottom=0.11,left=0.15,right=0.9)
            plt.text(x=1.5, y=1.5, verticalalignment='center', s='Source', size=12, rotation=-90)        
            plt.text(x=1.5, y=4.5, verticalalignment='center', s='Pivot', size=12, rotation=-90)        
            plt.text(x=1.5, y=7.5, verticalalignment='center', s='Sink', size=12, rotation=-90)        
        if sf == 'NEE':
            plt.subplots_adjust(top=0.95,bottom=0.11,left=0.33,right=0.99)
            plt.text(x=-1.45, y=1.3, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase', size=12, rotation=90)        
            
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_P1_bias_variance_phase_sub_box_whisker.png')
        plt.close()
        
        
#        Previous plot with 2 source sites and box plots
#        sns.set(font_scale=1.25, style='white')
#        g = sns.catplot(x='optim', y='value', hue="variable", row='site_type',
#                data=msd_melt_p1, kind="box", size=2, height=1, aspect=1, margin_titles=True, legend=False)
#        g.set(xlabel='', ylabel='')
#        g = g.map(plt.axvline, x=0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
#        [plt.setp(ax.texts, text="") for ax in g.axes.flat]
#                
#        plt.subplots_adjust(top=0.93,bottom=0.1,left=0.35,right=0.99)
#        if sf == 'Reco':
#            plt.legend(title='', bbox_to_anchor=(0.15,1.3), loc=2, fontsize=11)
#            g.set_titles(row_template = '{row_name}')
#            plt.subplots_adjust(top=0.93,bottom=0.1,left=0.16,right=0.87)
#        
#        plt.suptitle(subplot_title[nf]+' '+sf, fontsize=14)
#        
#        if sf == 'NEE':
#            plt.text(x=-1.5, y=3.5*max(msd_melt_p1.value)/2, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase', size=14, rotation=90)        
#        g.set_xticklabels(df_param_group.param[[0,8]], rotation=90, fontsize=14)
#        
#        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_P1_bias_variance_phase_sub_box_whisker.png')
#        plt.close()
        
        
        ## ---------------------
        #Box-whisker plots for mean annual flux difference
        df_annual_all_plot = df_annual_all.loc[df_annual_all['Flux'] == sf]
        
        df_annual_all_plot['Obs_Prior'] = df_annual_all_plot.Obs - df_annual_all_plot.Prior
        df_annual_all_plot['Obs_Post'] = df_annual_all_plot.Obs - df_annual_all_plot.Post
        
        flux_diff = df_annual_all_plot.groupby(['Site','Flux','optim','optim_name'])[['Obs_Post']].mean()
        flux_diff['Obs_Prior'] = df_annual_all_plot.groupby(['Site','Flux'])[['Obs_Prior']].mean()
        flux_diff = flux_diff.reset_index()
        
        df_flux_diff = pd.melt(flux_diff, id_vars=['Site','Flux','optim','optim_name'], var_name="Sim", value_name="Flux_diff")
        df_flux_diff.optim[df_flux_diff.Sim == 'Obs_Prior'] = 0
        df_flux_diff.optim_name[df_flux_diff.Sim == 'Obs_Prior'] = 'Prior'
        df_flux_diff.drop_duplicates(keep='first',inplace=True) # dropping duplicate values 
        
        g = sns.boxplot(x='optim', y='Flux_diff', hue='optim_name', data=df_flux_diff, 
                        order=boxplot_order, hue_order=['Prior','nee','gpp_reco','gpp','reco'], dodge=False)  
        if sf == 'NEE':
            plt.axhline(y=0, ls='--', c='grey') # add a horizontal line at ratio = 1.0
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
        plt.text(x=-1.6, y=(max(df_flux_diff.Flux_diff)+min(df_flux_diff.Flux_diff))/2, verticalalignment='center', s='Flux mean annual deviation ($gCm^{-2}year^{-1}$)', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param, rotation=90, fontsize=10)
        new_labels = ['Prior','NEE','GPP+$R_{eco}$','GPP','$R_{eco}$']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_diff_box_whisker.png')
        plt.close()
        
        #Box-whisker plots for mean annual flux difference (only NEE opt)
        df_flux_diff = df_flux_diff[df_flux_diff['optim_name'].isin(['Prior','nee'])]
        
        sns.set(font_scale=1.1, style='white')
        g = sns.boxplot(x='optim', y='Flux_diff', hue='optim_name', data=df_flux_diff, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        if sf == 'NEE':
            plt.axhline(y=0, ls='--', c='grey') # add a horizontal line at ratio = 1.0
        g.figure.set_size_inches(5,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.99,bottom=0.12,left=0.15,right=0.99)
        plt.legend(title='Optimization', loc='best', fontsize=14)
#        plt.suptitle(sf+' Flux')   
        plt.text(x=3.3, y=170, s='(a)', size=16)  
        plt.text(x=-1.9, y=(max(df_flux_diff.Flux_diff)+min(df_flux_diff.Flux_diff))/2, verticalalignment='center', s='Flux mean annual deviation ($gCm^{-2}year^{-1}$)', size=14, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_diff_sub_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for bias, variance and phase fractions side by side (only NEE opt)
        flux_msd_all_plot = flux_msd_all.loc[flux_msd_all['Flux'] == sf]
        
        flux_msd_all_plot.optim[flux_msd_all_plot.Sim == 'prior'] = 0
        flux_msd_all_plot.optim_name[flux_msd_all_plot.Sim == 'prior'] = 'Prior'
        
        flux_msd_all_plot = flux_msd_all_plot.drop('r', 1)
        flux_msd_all_plot = pd.merge(flux_msd_all_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
        
        #calculate the fractions in terms of total MSD        
        flux_msd_all_plot.bias = flux_msd_all_plot.bias / flux_msd_all_plot.mse
        flux_msd_all_plot.variance = flux_msd_all_plot.variance / flux_msd_all_plot.mse
        flux_msd_all_plot.phase = flux_msd_all_plot.phase / flux_msd_all_plot.mse
        
        flux_msd_all_plot = flux_msd_all_plot[flux_msd_all_plot['optim_name'].isin(['Prior','nee'])]
        msd_melt = pd.melt(flux_msd_all_plot, id_vars=('Sim','Site','site_type','Flux','optim','optim_name','param'), var_name="variable", value_name="value")
        
        g = sns.boxplot(x='optim', y='value', hue='variable', data=msd_melt, order=boxplot_order_sub,
                      hue_order=['bias','variance','phase'])  
        g.figure.set_size_inches(4,4)
        g.set(xlabel='', ylabel='')
        plt.axvline(x=0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g.legend_.remove()
                
        plt.subplots_adjust(top=0.93,bottom=0.12,left=0.16,right=0.99)
        if sf == 'Reco':
            plt.legend(title='', loc='upper center', fontsize=11)
            
        plt.suptitle(sf+' Flux')   
        if sf == 'NEE':
            plt.text(x=-1.9, y=0.5, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase contribution to MSD', size=12, rotation=90)        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=11)
        
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_daily_bias_variance_phase_sub_fraction_box_whisker.png')
        plt.close()
        
        
    #%% plotting MSD comparison for monthly time series
    ## ---------------------
    ## ---------------------
    #assign site type (source / pivot / sink)
    flux_msd_all_month['site_type'] = 'Sink'
#    flux_msd_all.loc[flux_msd_all['Site'] == 'US-SO2', 'site_type'] = 'Pivot'
    flux_msd_all_month.loc[flux_msd_all_month['Site'] == 'US-Wkg', 'site_type'] = 'Pivot'
    flux_msd_all_month.loc[flux_msd_all_month['Site'] == 'US-Whs', 'site_type'] = 'Pivot'
    flux_msd_all_month.loc[flux_msd_all_month['Site'] == 'US-SRG', 'site_type'] = 'Pivot'
    flux_msd_all_month.loc[flux_msd_all_month['Site'] == 'US-SRM', 'site_type'] = 'Pivot'
    flux_msd_all_month.loc[flux_msd_all_month['Site'] == 'US-Seg', 'site_type'] = 'Pivot'
        
#    flux_msd_all.loc[flux_msd_all['Site'] == 'US-Var', 'site_type'] = 'Source'
#    flux_msd_all_month.loc[flux_msd_all_month['Site'] == 'US-Seg', 'site_type'] = 'Source'
    flux_msd_all_month.loc[flux_msd_all_month['Site'] == 'US-Aud', 'site_type'] = 'Source'

    for nf, sf in enumerate(flux_names):
#        nf = 0
        sf = flux_names[nf]
        
        flux_msd_all_month_plot = flux_msd_all_month.loc[flux_msd_all_month['Flux'] == sf]
        
        flux_msd_all_month_plot.optim[flux_msd_all_month_plot.Sim == 'prior'] = 0
        flux_msd_all_month_plot.optim_name[flux_msd_all_month_plot.Sim == 'prior'] = 'Prior'
        
        flux_msd_all_month_plot = flux_msd_all_month_plot.drop('mse', 1)
        flux_msd_all_month_plot = flux_msd_all_month_plot.drop('r', 1)
        flux_msd_all_month_plot = pd.merge(flux_msd_all_month_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
        
        
        ## ---------------------
        #Box-whisker plots for bias (only NEE opt)
        flux_msd_all_month_plot = flux_msd_all_month_plot[flux_msd_all_month_plot['optim_name'].isin(['Prior','nee'])]
        
        g = sns.boxplot(x='optim', y='bias', hue='optim_name', data=flux_msd_all_month_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
        plt.text(x=-0.95, y=max(flux_msd_all_plot.bias)/2, verticalalignment='center', s='$Bias^{2}$', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_monthly_bias_sub_box_whisker.png')
        plt.close()
        
        #Box-whisker plots for variance (only NEE opt)
        g = sns.boxplot(x='optim', y='variance', hue='optim_name', data=flux_msd_all_month_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
        plt.text(x=-0.95, y=max(flux_msd_all_plot.variance)/2, verticalalignment='center', s='$Variance^{2}$', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_monthly_variance_sub_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for phase (only NEE opt)
        g = sns.boxplot(x='optim', y='phase', hue='optim_name', data=flux_msd_all_month_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', bbox_to_anchor=(0.6,1))
        plt.suptitle(sf+' Flux')   
        plt.text(x=-0.95, y=max(flux_msd_all_plot.phase)/2, verticalalignment='center', s='Phase', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_monthly_phase_sub_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for bias, variance and phase side by side (only NEE opt)
        msd_month_melt = pd.melt(flux_msd_all_month_plot, id_vars=('Sim','Site','site_type','Flux','optim','optim_name','param'), var_name="variable", value_name="value")
        
        sns.set(font_scale=1.25, style='white')
        g = sns.catplot(x='optim', y='value', hue="variable", row='site_type',
                data=msd_month_melt, kind="box", size=3, height=1, aspect=1.5, margin_titles=True, legend=False, sharex=False, sharey=False)
        g.set(xlabel='', ylabel='')
        g = g.map(plt.axvline, x=0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        [plt.setp(ax.texts, text="") for ax in g.axes.flat]
                
        plt.subplots_adjust(top=0.95,bottom=0.06,left=0.16,right=0.99)
        if sf == 'Reco':
            plt.text(x=7, y=11.4, s='(c)', size=20)  
            plt.legend(title='', loc='upper center', fontsize=11)
            g.set_titles(row_template = '{row_name}')
            plt.subplots_adjust(top=0.95,bottom=0.07,left=0.16,right=0.93)
        
        plt.suptitle(sf+' Flux')
        if sf == 'GPP':
            plt.text(x=6.5, y=7.25, s='(b)', size=20)  
        
        if sf == 'NEE':
            plt.text(x=6.5, y=3.41, s='(a)', size=20)  
            plt.text(x=-1.9, y=3.5*max(msd_month_melt.value)/2, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase', size=14, rotation=90)        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
        
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_monthly_bias_variance_phase_sub_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for bias, variance and phase fractions side by side (only NEE opt)
        flux_msd_all_month_plot = flux_msd_all_month.loc[flux_msd_all_month['Flux'] == sf]
        
        flux_msd_all_month_plot.optim[flux_msd_all_month_plot.Sim == 'prior'] = 0
        flux_msd_all_month_plot.optim_name[flux_msd_all_month_plot.Sim == 'prior'] = 'Prior'
        
        flux_msd_all_month_plot = flux_msd_all_month_plot.drop('r', 1)
        flux_msd_all_month_plot = pd.merge(flux_msd_all_month_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
        
        #calculate the fractions in terms of total MSD        
        flux_msd_all_month_plot.bias = flux_msd_all_month_plot.bias / flux_msd_all_month_plot.mse
        flux_msd_all_month_plot.variance = flux_msd_all_month_plot.variance / flux_msd_all_month_plot.mse
        flux_msd_all_month_plot.phase = flux_msd_all_month_plot.phase / flux_msd_all_month_plot.mse
        
        flux_msd_all_month_plot = flux_msd_all_month_plot[flux_msd_all_month_plot['optim_name'].isin(['Prior','nee'])]
        msd_month_melt = pd.melt(flux_msd_all_month_plot, id_vars=('Sim','Site','site_type','Flux','optim','optim_name','param'), var_name="variable", value_name="value")
        
        g = sns.boxplot(x='optim', y='value', hue='variable', data=msd_month_melt, order=boxplot_order_sub,
                      hue_order=['bias','variance','phase'])  
        g.figure.set_size_inches(4,4)
        g.set(xlabel='', ylabel='')
        plt.axvline(x=0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g.legend_.remove()
                
        plt.subplots_adjust(top=0.93,bottom=0.12,left=0.16,right=0.99)
        if sf == 'Reco':
            plt.legend(title='', loc='upper center', fontsize=11)
            
        plt.suptitle(sf+' Flux')   
        if sf == 'NEE':
            plt.text(x=-1.9, y=0.5, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase contribution to MSD', size=12, rotation=90)        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=11)
        
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_monthly_bias_variance_phase_sub_fraction_box_whisker.png')
        plt.close()
        
        
    #%% plotting MSD comparison for annual time series
    ## ---------------------
    ## ---------------------
    #assign site type (source / pivot / sink)
    flux_msd_all_annual['site_type'] = 'Sink'
#    flux_msd_all.loc[flux_msd_all['Site'] == 'US-SO2', 'site_type'] = 'Pivot'
    flux_msd_all_annual.loc[flux_msd_all_annual['Site'] == 'US-Wkg', 'site_type'] = 'Pivot'
    flux_msd_all_annual.loc[flux_msd_all_annual['Site'] == 'US-Whs', 'site_type'] = 'Pivot'
    flux_msd_all_annual.loc[flux_msd_all_annual['Site'] == 'US-SRG', 'site_type'] = 'Pivot'
    flux_msd_all_annual.loc[flux_msd_all_annual['Site'] == 'US-SRM', 'site_type'] = 'Pivot'
    flux_msd_all_annual.loc[flux_msd_all_annual['Site'] == 'US-Seg', 'site_type'] = 'Pivot'
        
#    flux_msd_all.loc[flux_msd_all['Site'] == 'US-Var', 'site_type'] = 'Source'
    flux_msd_all_annual.loc[flux_msd_all_annual['Site'] == 'US-Aud', 'site_type'] = 'Source'

    for nf, sf in enumerate(flux_names):
#        nf = 0
        sf = flux_names[nf]
        
        flux_msd_all_annual_plot = flux_msd_all_annual.loc[flux_msd_all_annual['Flux'] == sf]
        
        flux_msd_all_annual_plot.optim[flux_msd_all_annual_plot.Sim == 'prior'] = 0
        flux_msd_all_annual_plot.optim_name[flux_msd_all_annual_plot.Sim == 'prior'] = 'Prior'
        
        flux_msd_all_annual_plot = flux_msd_all_annual_plot.drop('mse', 1)
        flux_msd_all_annual_plot = flux_msd_all_annual_plot.drop('r', 1)
        flux_msd_all_annual_plot = pd.merge(flux_msd_all_annual_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
        
        
        ## ---------------------
        #Box-whisker plots for bias (only NEE opt)
        flux_msd_all_annual_plot = flux_msd_all_annual_plot[flux_msd_all_annual_plot['optim_name'].isin(['Prior','nee'])]
        
        g = sns.boxplot(x='optim', y='bias', hue='optim_name', data=flux_msd_all_annual_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
        plt.text(x=-0.95, y=max(flux_msd_all_plot.bias)/2, verticalalignment='center', s='$Bias^{2}$', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_bias_sub_box_whisker.png')
        plt.close()
        
        #Box-whisker plots for variance (only NEE opt)
        g = sns.boxplot(x='optim', y='variance', hue='optim_name', data=flux_msd_all_annual_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', loc='best')
        plt.suptitle(sf+' Flux')   
        plt.text(x=-0.95, y=max(flux_msd_all_plot.variance)/2, verticalalignment='center', s='$Variance^{2}$', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_variance_sub_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for phase (only NEE opt)
        g = sns.boxplot(x='optim', y='phase', hue='optim_name', data=flux_msd_all_annual_plot, 
                        order=boxplot_order_sub, hue_order=['Prior','nee'], dodge=False)  
        g.figure.set_size_inches(10,5)
        g.set(xlabel='', ylabel='')
        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.07,right=0.99)
        plt.legend(title='Optimization', bbox_to_anchor=(0.6,1))
        plt.suptitle(sf+' Flux')   
        plt.text(x=-0.95, y=max(flux_msd_all_plot.phase)/2, verticalalignment='center', s='Phase', size=12, rotation=90)
        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=10)
        new_labels = ['Prior','NEE']
        leg = g.axes.get_legend()
        for t, l in zip(leg.texts, new_labels): t.set_text(l)

        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_phase_sub_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for bias, variance and phase side by side (only NEE opt)
        msd_annual_melt = pd.melt(flux_msd_all_annual_plot, id_vars=('Sim','Site','site_type','Flux','optim','optim_name','param'), var_name="variable", value_name="value")
        msd_annual_melt = msd_annual_melt.drop_duplicates(subset=['Sim','Site','site_type','Flux','optim','optim_name','param','variable'], keep="first")

        sns.set(font_scale=1.25, style='white')
        g = sns.catplot(x='optim', y='value', hue="variable", row='site_type',
                data=msd_annual_melt, kind="box", size=3, height=1, aspect=1.5, margin_titles=True, legend=False, sharex=False, sharey=False)
        g.set(xlabel='', ylabel='')
        g = g.map(plt.axvline, x=0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        [plt.setp(ax.texts, text="") for ax in g.axes.flat]
                
        plt.subplots_adjust(top=0.95,bottom=0.07,left=0.21,right=0.99)
        if sf == 'Reco':
            plt.text(x=7, y=11.4, s='(c)', size=20)  
            plt.legend(title='', loc='upper center', fontsize=11)
            g.set_titles(row_template = '{row_name}')
            plt.subplots_adjust(top=0.95,bottom=0.06,left=0.16,right=0.93)
        
        plt.suptitle(sf+' Flux')
        if sf == 'GPP':
            plt.text(x=6.5, y=7.25, s='(b)', size=20)  
        
        if sf == 'NEE':
            plt.text(x=6.5, y=3.41, s='(a)', size=20)  
            plt.text(x=-2.5, y=3.5*max(msd_annual_melt.value)/2, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase', size=14, rotation=90)        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=14)
        
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_bias_variance_phase_sub_box_whisker.png')
        plt.close()
        
        #Box-whisker plots for bias, variance and phase side by side (only NEE opt and Prior + P1)
        msd_annual_melt_p1 = msd_annual_melt[msd_annual_melt['param'].isin(['Prior','P1'])]
        
        sns.set(font_scale=1, style='white')
        fig, axes = plt.subplots(3, 1, figsize=(2,5))
        axes = axes.flatten()
        g1 = sns.boxplot(x='optim', y='value', hue='variable', data=msd_annual_melt_p1[msd_annual_melt_p1['site_type'].isin(['Sink'])], dodge=True, ax=axes[0])  
        g1.set(xlabel='', ylabel='')
        g1.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g1.set(xticklabels=[])
        g1.legend_.remove()
#        g1.set(title = subplot_title[nf]+' '+sf)
        g1.set(title = sf)
        if sf == 'NEE':
            g1.set(ylim=(-0.03*max(msd_annual_melt_p1.value), 0.03*max(msd_annual_melt_p1.value)+max(msd_annual_melt_p1.value)))
            g1.text(x=1.1, y=85000, verticalalignment='center', s='(a)', size=12)        
#        else: g1.set(ylim=(-10000,390000))
        if sf == 'GPP': g1.text(x=1.1, y=170000, verticalalignment='center', s='(b)', size=12)        
        if sf == 'Reco': g1.text(x=1.1, y=360000, verticalalignment='center', s='(c)', size=12)        
                
        g2 = sns.boxplot(x='optim', y='value', hue='variable', data=msd_annual_melt_p1[msd_annual_melt_p1['site_type'].isin(['Pivot'])], dodge=True, ax=axes[1])  
        g2.set(xlabel='', ylabel='')
        g2.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g2.set(xticklabels=[])
        g2.legend_.remove()
        if sf == 'NEE':
#            g2.set(ylim=(-0.03*max(msd_annual_melt_p1.value), 0.03*max(msd_annual_melt_p1.value)+max(msd_annual_melt_p1.value)))
            g2.set(ylim=(-50, 2800))
            g2.text(x=1.1, y=2400, verticalalignment='center', s='(d)', size=12)        
#        else: g2.set(ylim=(-10000,390000))
        if sf == 'GPP': g2.text(x=1.1, y=21000, verticalalignment='center', s='(e)', size=12)        
        if sf == 'Reco': g2.text(x=1.1, y=19000, verticalalignment='center', s='(f)', size=12)        
        
        g3 = sns.barplot(x="optim", y="value", hue="variable", data=msd_annual_melt_p1[msd_annual_melt_p1['site_type'].isin(['Source'])], dodge=True, ax=axes[2])
#        g3 = sns.stripplot(x="optim", y="value", edgecolor="none", hue="variable", linewidth=2,
#                  marker='_', data=msd_melt_p1[msd_melt_p1['site_type'].isin(['Source'])], dodge=True, ax=axes[2])
        g3.set(xlabel='', ylabel='')
        g3.axvline(0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g3.set_xticklabels(df_param_group.param[0:2], rotation=90, fontsize=12)
        g3.legend_.remove()
        if sf == 'NEE':
            g3.set(ylim=(-0.03*max(msd_annual_melt_p1.value), 0.03*max(msd_annual_melt_p1.value)+max(msd_annual_melt_p1.value)))
            g3.set(ylim=(-500, 30000))
            g3.text(x=1.1, y=27000, verticalalignment='center', s='(g)', size=12)        
#        else: g3.set(ylim=(-10000,390000))
        if sf == 'GPP': g3.text(x=1.1, y=10500, verticalalignment='center', s='(h)', size=12)        
        if sf == 'Reco': g3.text(x=1.1, y=6000, verticalalignment='center', s='(i)', size=12)        
        
        plt.subplots_adjust(wspace=0, hspace=0.03)
        if sf == 'GPP': plt.subplots_adjust(top=0.95,bottom=0.11,left=0.33,right=0.99)
        if sf == 'Reco':
            plt.legend(title='', loc='best', fontsize=10)
            plt.subplots_adjust(top=0.95,bottom=0.11,left=0.33,right=0.9)
            plt.text(x=1.5, y=3000, verticalalignment='center', s='Source', size=12, rotation=-90)        
            plt.text(x=1.5, y=10000, verticalalignment='center', s='Pivot', size=12, rotation=-90)        
            plt.text(x=1.5, y=18000, verticalalignment='center', s='Sink', size=12, rotation=-90)        
        if sf == 'NEE':
            plt.subplots_adjust(top=0.95,bottom=0.11,left=0.38,right=0.99)
            plt.text(x=-1.7, y=50000, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase', size=12, rotation=90)        
            
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_P1_bias_variance_phase_sub_box_whisker.png')
        plt.close()
        
        
        #Box-whisker plots for bias, variance and phase fractions side by side (only NEE opt)
        flux_msd_all_annual_plot = flux_msd_all_annual.loc[flux_msd_all_annual['Flux'] == sf]
        
        flux_msd_all_annual_plot.optim[flux_msd_all_annual_plot.Sim == 'prior'] = 0
        flux_msd_all_annual_plot.optim_name[flux_msd_all_annual_plot.Sim == 'prior'] = 'Prior'
        
        flux_msd_all_annual_plot = flux_msd_all_annual_plot.drop('r', 1)
        flux_msd_all_annual_plot = pd.merge(flux_msd_all_annual_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
        
        #calculate the fractions in terms of total MSD        
        flux_msd_all_annual_plot.bias = flux_msd_all_annual_plot.bias / flux_msd_all_annual_plot.mse
        flux_msd_all_annual_plot.variance = flux_msd_all_annual_plot.variance / flux_msd_all_annual_plot.mse
        flux_msd_all_annual_plot.phase = flux_msd_all_annual_plot.phase / flux_msd_all_annual_plot.mse
        
        flux_msd_all_annual_plot = flux_msd_all_annual_plot[flux_msd_all_annual_plot['optim_name'].isin(['Prior','nee'])]
        msd_annual_melt = pd.melt(flux_msd_all_annual_plot, id_vars=('Sim','Site','site_type','Flux','optim','optim_name','param'), var_name="variable", value_name="value")
        
        g = sns.boxplot(x='optim', y='value', hue='variable', data=msd_annual_melt, order=boxplot_order_sub,
                      hue_order=['bias','variance','phase'])  
        g.figure.set_size_inches(4,4)
        g.set(xlabel='', ylabel='')
        plt.axvline(x=0.5, ls='--', c='black') # add a horizontal line at ratio = 1.0
        g.legend_.remove()
                
        plt.subplots_adjust(top=0.93,bottom=0.12,left=0.16,right=0.99)
        if sf == 'Reco':
            plt.legend(title='', loc='upper center', fontsize=11)
            
        plt.suptitle(sf+' Flux')   
        if sf == 'NEE':
            plt.text(x=-1.9, y=0.5, verticalalignment='center', s='$Bias^{2}$ / $Variance^{2}$ / Phase contribution to MSD', size=12, rotation=90)        
        g.set_xticklabels(df_param_group.param[0:8], rotation=90, fontsize=11)
        
        plt.savefig(outdir + 'Flux_annual_anomaly/box_whisker_plots/' + sf + '_flux_annual_bias_variance_phase_sub_fraction_box_whisker.png')
        plt.close()
        
        
#%% Scatter plots for annual fluxes
if scatter_plot == True: 
    optim1 = ['(GPP + $R_{eco}$)','NEE','GPP','$R_{eco}$']
    for nf, sf in enumerate(flux_names):
#        nf = 0
        sf = flux_names[nf]
        
        for nj, nf in enumerate(df_annual_all['optim_name'].unique()):
#            nj = 1
            nf = df_annual_all['optim_name'].unique()[nj]
        
            df_annual_all_plot = df_annual_all.drop('Date', 1)
            df_annual_all_plot = pd.melt(df_annual_all_plot, id_vars=('Year','Obs','Site','Flux','optim','optim_name'), var_name="Simulation", value_name="value")
            df_annual_all_plot = df_annual_all_plot.loc[df_annual_all_plot['Flux'] == sf]
            
            df_annual_all_plot.optim[df_annual_all_plot.Simulation == 'Prior'] = 0
            df_annual_all_plot.optim_name[df_annual_all_plot.Simulation == 'Prior'] = 'Prior'
            
            df_annual_all_plot = pd.merge(df_annual_all_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
            
            df_annual_all_plot = df_annual_all_plot[df_annual_all_plot['optim_name'].isin(['Prior',nf])]
            df_annual_all_plot = df_annual_all_plot.drop_duplicates(subset=['Year','Site','Flux','optim','optim_name','Simulation','param'], keep="first")
        
#            df_annual_all_sub = df_annual_all_plot[df_annual_all_plot['optim_name'].isin([nf])]
#            df_annual_all_sub = df_annual_all_sub[df_annual_all_sub['Site'].isin(['US-Aud'])]
#            df_annual_all_sub = df_annual_all_sub.drop_duplicates(subset=['Year','Site','Flux','optim','optim_name','Simulation','param'], keep="first")
        
            # fit a model with all interactions
            fit = smf.ols('value ~ Obs * param * Site', df_annual_all_plot).fit()
            df_annual_all_plot['yhat'] = fit.predict(df_annual_all_plot[['Obs','param','Site']])
    
            sns.set(font_scale=1.6, style='white')
            colors = sns.color_palette("hls", 13)
            g = sns.FacetGrid(df_annual_all_plot, col="param", hue="Site", margin_titles=True,
                           col_wrap = 4, sharex=False, sharey=False, height=12, aspect=0.5, palette=colors)
            g.map(plt.plot, "Obs", "yhat", linewidth=2)
            #            g.add_legend(bbox_to_anchor=(0.16,0.86), ncol=2, fontsize=11)
            g.add_legend(bbox_to_anchor=(1,0.5), ncol=1, fontsize=14) 
            g.map(plt.scatter, "Obs", "value", marker='.', s=35).set(xlim=(-450,450) , ylim=(-450,450))

            axes = g.fig.axes
            x = np.arange(-450, 450, 50)
            y = x
            for ax in axes:
                ax.plot(y, x, C='grey', linestyle="--")
                ax.axhline(0, C='k', ls='--')
                ax.axvline(0, C='k', ls='--')
            
            g.fig.set_size_inches(18,10)
            g.set(xlabel='', ylabel='')
            [plt.setp(ax.texts, text="") for ax in g.axes.flat]
            g.set_titles(col_template = '{col_name}')
        
            plt.subplots_adjust(top=0.96,bottom=0.09,left=0.07,right=0.9)
#            plt.suptitle(sf+' Flux')   
            plt.text(x=-4230, y=600, verticalalignment='center', s='Model NEE ($gCm^{-2}year^{-1}$)', size=20, rotation=90)
            plt.text(x=-1800, y=-600, verticalalignment='center', s='Site NEE ($gCm^{-2}year^{-1}$)', size=20, rotation=0)
            
            plt.savefig(outdir + 'Scatter_plots/' + sf + '_scatter_plots_for_' + optim1[nj] + '_optimization.png')
            plt.close()
     
        
    #Correlation coefficient R
    flux_corr_all_plot = flux_corr_all.loc[flux_corr_all['Flux'] == 'NEE']
    flux_corr_all_plot.optim[flux_corr_all_plot.Simulation == 'Prior'] = 0
    flux_corr_all_plot.optim_name[flux_corr_all_plot.Simulation == 'Prior'] = 'Prior'
    flux_corr_all_plot = pd.merge(flux_corr_all_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
    flux_corr_all_plot = flux_corr_all_plot.loc[flux_corr_all_plot['optim'].isin([0,8])]
    flux_corr_all_plot = flux_corr_all_plot.drop_duplicates(['Site','Simulation','Flux','optim','optim_name','param'], keep='first')
        
    #Slope of regression
    flux_slope_all_plot = flux_slope_all.loc[flux_slope_all['Flux'] == 'NEE']
    flux_slope_all_plot.optim[flux_slope_all_plot.Simulation == 'Prior'] = 0
    flux_slope_all_plot.optim_name[flux_slope_all_plot.Simulation == 'Prior'] = 'Prior'   
    flux_slope_all_plot = pd.merge(flux_slope_all_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
    flux_slope_all_plot = flux_slope_all_plot.loc[flux_slope_all_plot['optim'].isin([0,8])]
    flux_slope_all_plot = flux_slope_all_plot.drop_duplicates(['Site','Simulation','Flux','optim','optim_name','param'], keep='first')
    
    
    ## Scatter plots for annual fluxes for only prior and Opt1 P1
    df_annual_all_plot = df_annual_all.drop('Date', 1)
    df_annual_all_plot = pd.melt(df_annual_all_plot, id_vars=('Year','Obs','Site','Flux','optim','optim_name'), var_name="Simulation", value_name="value")
    df_annual_all_plot = df_annual_all_plot.loc[df_annual_all_plot['Flux'] == 'NEE']
            
    df_annual_all_plot.optim[df_annual_all_plot.Simulation == 'Prior'] = 0
    df_annual_all_plot.optim_name[df_annual_all_plot.Simulation == 'Prior'] = 'Prior'
            
    df_annual_all_plot = pd.merge(df_annual_all_plot,df_param_group,left_on='optim',right_on='optim',how='outer')
            
#    df_annual_all_plot = df_annual_all_plot[df_annual_all_plot['optim_name'].isin(['Prior','nee'])]
    df_annual_all_plot = df_annual_all_plot[df_annual_all_plot['optim'].isin([0,8])]
    df_annual_all_plot = df_annual_all_plot.drop_duplicates(subset=['Year','Site','Flux','optim','optim_name','Simulation','param'], keep="first")
            
    # fit a model with all interactions
    fit = smf.ols('value ~ Obs * param * Site', df_annual_all_plot).fit()
    df_annual_all_plot['yhat'] = fit.predict(df_annual_all_plot[['Obs','param','Site']])
    
    sns.set(font_scale=1.2, style='white')
    colors = sns.color_palette("hls", 13)
    g = sns.FacetGrid(df_annual_all_plot, col="param", hue="Site", margin_titles=True,
        sharex=False, sharey=False, height=12, aspect=0.5, palette=colors)
    g.map(plt.plot, "Obs", "yhat", linewidth=2)
    g.add_legend(bbox_to_anchor=(1,0.5), ncol=1, fontsize=14)
    
    g.map(plt.scatter, "Obs", "value", marker='.', s=35).set(xlim=(-450,450) , ylim=(-450,450))
#    g.add_legend(bbox_to_anchor=(0.27,0.8), ncol=2, fontsize=10)
    g.fig.set_size_inches(13,6)
            
    axes = g.fig.axes
    x = np.arange(-450, 500, 50)
    y = x
    for ax in axes:
        ax.plot(y, x, C='grey', linestyle="--")
        ax.axhline(0, C='k', ls='--')
        ax.axvline(0, C='k', ls='--')
            
    for ax in g.axes.flatten(): # Loop directly on the flattened axes 
        for _, spine in ax.spines.items():
            spine.set_visible(True) # You have to first turn them on
            spine.set_color('black')
            spine.set_linewidth(2)
        
    g.set(xlabel='', ylabel='')
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    axes = g.axes.flatten()
    axes[0].set_title("")
    axes[1].set_title("")
    plt.text(x=-990, y=400, s='(a) Prior', size=20)
    plt.text(x=30, y=400, s='(b) Posterior', size=20)
    
    plt.subplots_adjust(top=0.99,bottom=0.11,left=0.07,right=0.88)
    plt.suptitle('')   
    plt.text(x=-1620, y=0, verticalalignment='center', s='Model NEE ($gCm^{-2}year^{-1}$)', size=18, rotation=90)
    plt.text(x=-680, y=-530, verticalalignment='center', s='Site NEE ($gCm^{-2}year^{-1}$)', size=18, rotation=0)
    
#    df = pd.DataFrame(np.random.randint(low=0, high=10, size=(5, 3)),
#                    columns=['a', 'b', 'c'])
    
#    ycord = np.linspace(420,70,12)
#    par = flux_slope_all_plot.param.unique()[0]
#    flux_slope_all_plot_sub = flux_slope_all_plot.loc[flux_slope_all_plot['param'] == par]
#    flux_corr_all_plot_sub = flux_corr_all_plot.loc[flux_corr_all_plot['param'] == par]
#    for k, val in enumerate(flux_slope_all_plot_sub['Anomaly']):
#        r = str(round(flux_corr_all_plot_sub.Anomaly.iloc[k], 3))
#        m = str(round(flux_slope_all_plot_sub.Anomaly.iloc[k], 3))
#        plt.text(-1440, ycord[k], s='r = '+str(r)+'; m = '+str(m), size=12, color=colors[k])
#    
#    par = flux_slope_all_plot.param.unique()[1]
#    flux_slope_all_plot_sub = flux_slope_all_plot.loc[flux_slope_all_plot['param'] == par]
#    flux_corr_all_plot_sub = flux_corr_all_plot.loc[flux_corr_all_plot['param'] == par]
#    for k, val in enumerate(flux_slope_all_plot_sub['Anomaly']):
#        r = str(round(flux_corr_all_plot_sub.Anomaly.iloc[k], 3))
#        m = str(round(flux_slope_all_plot_sub.Anomaly.iloc[k], 3))
#        plt.text(-420, ycord[k], s='r = '+str(r)+'; m = '+str(m), size=12, color=colors[k])
        
    plt.savefig(outdir + 'Scatter_plots/' + 'NEE_prior_P1_scatter_plots_for_NEE_optimization.png')
    plt.close()
    
    
    #plot inset figures with mean NEE vs slopes of linear regression lines (both prior and P1)
    df_annual_obs = df_annual_all.loc[df_annual_all['optim'] == 1]
    df_annual_obs = df_annual_obs.loc[df_annual_obs['Flux'] == 'NEE']
    df_annual_obs = df_annual_obs[['Site','Obs']]
    df_annual_obs = df_annual_obs.groupby('Site', as_index=False)['Obs'].mean()
    
    
    flux_slope_all_plot_prior = flux_slope_all_plot.loc[flux_slope_all_plot['param'] == 'Prior']
    flux_slope_all_plot_prior.drop_duplicates(subset ="Site", inplace = True) 
    flux_corr_all_plot_prior = flux_corr_all_plot.loc[flux_corr_all_plot['param'] == 'Prior']
    flux_corr_all_plot_prior.drop_duplicates(subset ="Site", inplace = True) 
    
    flux_slope_all_plot_prior = pd.merge(flux_slope_all_plot_prior, df_annual_obs, how='left', on=['Site'])
    flux_corr_all_plot_prior = pd.merge(flux_corr_all_plot_prior, df_annual_obs, how='left', on=['Site'])
        
    # plot the inset figure for prior slope
    sns.set(font_scale=1.8, style='white')
    fig, ax = plt.subplots(figsize=(4,4))
    ax.bar(flux_slope_all_plot_prior.Site,flux_slope_all_plot_prior.Anomaly, color='grey')
#    ax.scatter(flux_slope_all_plot_prior.Obs,flux_slope_all_plot_prior.Anomaly)
#    plt.xticks([-300,-200,-100,0,100])
    fig.subplots_adjust(top=0.98,bottom=0.28,left=0.21,right=0.995)
    ax.axhline(1, C='k', ls='--')
    ax.set_ylabel('Slope')
#    ax.set_xlabel('Mean NEE ($gCm^{-2}year^{-1}$)')
    ax.set_ylim(-0.15, 1.25)
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.xticks(rotation=90, size=18)
    plt.savefig(outdir + 'Scatter_plots/' + 'NEE_prior_slope_inset_scatter_plots.png')
    plt.close()
    
    # plot the inset figure for prior slope
    fig, ax = plt.subplots(figsize=(4,4))
    ax.bar(flux_corr_all_plot_prior.Site,flux_corr_all_plot_prior.Anomaly, color='grey')
#    ax.scatter(flux_corr_all_plot_prior.Obs,flux_corr_all_plot_prior.Anomaly)
#    plt.xticks([-300,-200,-100,0,100])
    fig.subplots_adjust(top=0.99,bottom=0.28,left=0.26,right=0.995)
    ax.axhline(1, C='k', ls='--')
    ax.set_ylabel('Correlation (R)')
#    ax.set_xlabel('Mean NEE ($gCm^{-2}year^{-1}$)')
    ax.set_ylim(-0.5, 1.1)
    plt.yticks([-0.2,0,0.2,0.4,0.6,0.8,1])
    plt.xticks(rotation=90, size=18)
    plt.savefig(outdir + 'Scatter_plots/' + 'NEE_prior_corr_inset_scatter_plots.png')
    plt.close()
    
    flux_slope_all_plot_post = flux_slope_all_plot.loc[flux_slope_all_plot['param'] == 'P1']
    flux_corr_all_plot_post = flux_corr_all_plot.loc[flux_corr_all_plot['param'] == 'P1']
    flux_slope_all_plot_post = pd.merge(flux_slope_all_plot_post, df_annual_obs, how='left', on=['Site'])
    flux_corr_all_plot_post = pd.merge(flux_corr_all_plot_post, df_annual_obs, how='left', on=['Site'])
    
    # plot the inset figure for prior slope
    fig, ax = plt.subplots(figsize=(4,4))
    ax.bar(flux_slope_all_plot_post.Site,flux_slope_all_plot_post.Anomaly, color='grey')
#    ax.scatter(flux_slope_all_plot_post.Obs,flux_slope_all_plot_post.Anomaly)
#    plt.xticks([-300,-200,-100,0,100])
    fig.subplots_adjust(top=0.99,bottom=0.28,left=0.21,right=0.995)
    ax.axhline(1, C='k', ls='--')
    ax.set_ylabel('Slope')
#    ax.set_xlabel('Mean NEE ($gCm^{-2}year^{-1}$)')
    ax.set_ylim(-0.15, 1.25)
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.xticks(rotation=90, size=18)
    plt.savefig(outdir + 'Scatter_plots/' + 'NEE_P1_slope_inset_scatter_plots.png')
    plt.close()
    
    # plot the inset figure for prior slope
    fig, ax = plt.subplots(figsize=(4,4))
    ax.bar(flux_corr_all_plot_post.Site,flux_corr_all_plot_post.Anomaly, color='grey')
#    ax.scatter(flux_corr_all_plot_post.Obs,flux_corr_all_plot_post.Anomaly)
#    plt.xticks([-300,-200,-100,0,100])
    fig.subplots_adjust(top=0.99,bottom=0.28,left=0.26,right=0.995)
    ax.axhline(1, C='k', ls='--')
    ax.set_ylabel('Correlation (R)', size=20)
#    ax.set_xlabel('Mean NEE ($gCm^{-2}year^{-1}$)')
    ax.set_ylim(-0.5, 1.1)
    plt.yticks([-0.2,0,0.2,0.4,0.6,0.8,1])
    plt.xticks(rotation=90, size=18)
    plt.savefig(outdir + 'Scatter_plots/' + 'NEE_P1_corr_inset_scatter_plots.png')
    plt.close()
    
    
    
#%% ## Same scatter plots for NATASHA's NASA GRANT with annual fluxes for only prior and Opt1 P1
    sns.set(font_scale=1.2, style='white')
    colors = sns.color_palette("hls", 13)
    g = sns.FacetGrid(df_annual_all_plot, col="param", hue="Site", margin_titles=True,
        sharex=False, sharey=False, height=12, aspect=0.5, palette=colors)
    g.map(plt.plot, "Obs", "yhat")
    g.add_legend(bbox_to_anchor=(0.99,0.375), ncol=1, fontsize=12)
    
    g.map(plt.scatter, "Obs", "value", marker='.', s=25).set(xlim=(-450,450) , ylim=(-450,450))
#    g.add_legend(bbox_to_anchor=(0.27,0.8), ncol=2, fontsize=10)
    g.fig.set_size_inches(13,6)
            
    axes = g.fig.axes
    x = np.arange(-450, 500, 50)
    y = x
    for ax in axes:
        ax.plot(y, x, C='grey', linestyle="--")
        ax.axhline(0, C='k', ls='--')
        ax.axvline(0, C='k', ls='--')
            
    for ax in g.axes.flatten(): # Loop directly on the flattened axes 
        for _, spine in ax.spines.items():
            spine.set_visible(True) # You have to first turn them on
            spine.set_color('black')
            spine.set_linewidth(2)
        
    g.set(xlabel='', ylabel='')
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    axes = g.axes.flatten()
    axes[0].set_title("")
    axes[1].set_title("")
    plt.text(x=-990, y=400, s='(a) Prior', size=18)
    plt.text(x=30, y=400, s='(b) Posterior', size=18)
    
    plt.subplots_adjust(top=0.99,bottom=0.11,left=0.07,right=0.995)
    plt.suptitle('')   
    plt.text(x=-1600, y=0, verticalalignment='center', s='Model NEE ($gCm^{-2}year^{-1}$)', size=16, rotation=90)
    plt.text(x=-700, y=-530, verticalalignment='center', s='Site NEE ($gCm^{-2}year^{-1}$)', size=16, rotation=0)
    
    ycord = np.linspace(420,70,12)
    par = flux_slope_all_plot.param.unique()[0]
    flux_slope_all_plot_sub = flux_slope_all_plot.loc[flux_slope_all_plot['param'] == par]
    flux_corr_all_plot_sub = flux_corr_all_plot.loc[flux_corr_all_plot['param'] == par]
    for k, val in enumerate(flux_slope_all_plot_sub['Anomaly']):
#        r = str(round(flux_corr_all_plot_sub.Anomaly.iloc[k], 3))
        m = str(round(flux_slope_all_plot_sub.Anomaly.iloc[k], 3))
        plt.text(-1440, ycord[k], s='m = '+str(m), size=12, color=colors[k])
    
    par = flux_slope_all_plot.param.unique()[1]
    flux_slope_all_plot_sub = flux_slope_all_plot.loc[flux_slope_all_plot['param'] == par]
    flux_corr_all_plot_sub = flux_corr_all_plot.loc[flux_corr_all_plot['param'] == par]
    for k, val in enumerate(flux_slope_all_plot_sub['Anomaly']):
#        r = str(round(flux_corr_all_plot_sub.Anomaly.iloc[k], 3))
        m = str(round(flux_slope_all_plot_sub.Anomaly.iloc[k], 3))
        plt.text(-420, ycord[k], s='m = '+str(m), size=12, color=colors[k])
        
#    plt.legend(bbox_to_anchor=(0.5, 0.6), loc=1, borderaxespad=0.)
#    sns.plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
    
    plt.savefig(outdir + 'Scatter_plots/' + 'NEE_prior_P1_scatter_plots_for_NEE_optimization_NASA_Grant.png', dpi=150)
    plt.close()


#%% plotting RMSE comparison for gpp_reco vs nee optimization with all parameters
## ---------------------
## ---------------------
#make a table of all RMSE values for Opt1 P1 for Table S2
if rmse_plot == True: 
    df_rmse_all_plot = df_rmse_all.loc[df_rmse_all['Datatype'].isin(['Prior','Post'])]
    df_rmse_all_plot = df_rmse_all_plot.loc[df_rmse_all_plot['optim'].isin([8])]
    df_rmse_all_plot = df_rmse_all_plot.drop('optim', 1)
    df_rmse_all_plot = df_rmse_all_plot.drop('optim_name', 1)
    df_rmse_all_plot = df_rmse_all_plot.round(3) # round up the numbers to 3 decimal places
    
    df_rmse_pivot = df_rmse_all_plot.set_index(['Site','Datatype','Flux'])['RMSE'].unstack().reset_index()
    df_rmse_pivot = df_rmse_pivot.set_index(['Site','Datatype'])[(['GPP','Reco'])].unstack().reset_index()
    
    df_rmse_pivot = pd.merge(df_site_names,df_rmse_pivot,left_on='Site',right_on='Site',how='outer')
    
    #write the dataframe in outputs directory
    df_rmse_pivot.to_csv(outdir + 'RMSE/' + 'Flux_rmse_optim_' + optim_name[7] + '.csv', encoding='utf-8', index=False)
    


if rmse_plot == True and 8 in df_rmse_all.optim.unique(): 
    df_rmse_all_plot = df_rmse_all.loc[df_rmse_all['optim'].isin([1,8])]
    df_rmse_all_plot = df_rmse_all_plot.loc[df_rmse_all_plot['Datatype'] == 'Ratio']
    
    df_rmse_all_plot = df_rmse_all_plot.drop('Datatype', 1)
    df_rmse_all_plot = df_rmse_all_plot.drop('optim', 1)
              
    g = sns.catplot(x='optim_name', y='RMSE', hue='Flux', col='Site', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','blue','orange'], height=6, col_order = site_names, legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.89,left=0.03)
    plt.legend(loc='best')
    # entire figure title
    g.fig.suptitle('RMSE comparison '  + r"$\bf{WITH}$ " + param_title[0] + ': ' + optim[0] + ' optimization ' + r"$\bf{vs}$ " + optim[6] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
                                
    # save figure
    g.fig.savefig(outdir + 'RMSE/RMSE_comparison_gppreco_vs_nee_optimization.png')
    plt.close()


    ## plotting MEAN RMSE comparison for gpp_reco vs nee optimization with all parameters
    ## ---------------------       
    g = sns.catplot(x='optim_name', y='RMSE', hue='Flux', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','blue','orange'], legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.89,left=0.1,bottom=0.08)
    plt.legend(loc='best')
    # entire figure title
    g.fig.suptitle('Mean RMSE comparison '  + r"$\bf{WITH}$ " + param_title[0] + '\n' + optim[0] + ' optimization ' + r"$\bf{vs}$ " + optim[6] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
    g.fig.text(y=0.02, x=0.45, verticalalignment='center', s='Optimization', size=12, rotation=0)
    
    for ax in g.axes.flat:
        labels = ax.get_xticklabels() # get x labels
        labels[0] = optim[0] 
        labels[1] = optim[7]
        ax.set_xticklabels(labels) # set new labels 
                            
    # save figure
    g.fig.savefig(outdir + 'RMSE/Mean_RMSE_comparison_gppreco_vs_nee_optimization.png')
    plt.close()


    #%% plotting NEE optimization RMSE comparison for various parameter sets
    ## ---------------------
    ## ---------------------
if rmse_plot == True and [8,9,10,11,12,13,14] in df_rmse_all.optim.unique(): 
    df_rmse_all_plot = df_rmse_all.loc[df_rmse_all['optim'].isin([8,9,10,11,12,13,14])]
    df_rmse_all_plot = df_rmse_all_plot.loc[df_rmse_all_plot['Datatype'] == 'Ratio']
    
    df_rmse_all_plot = df_rmse_all_plot.drop('Datatype', 1)
    df_rmse_all_plot = df_rmse_all_plot.drop('optim_name', 1)
              
    g = sns.catplot(x='Flux', y='RMSE', hue='optim', col='Site', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','yellow','green','blue','orange','red'], height=6, 
                         col_wrap = col_wrap, col_order = site_names, legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.92,left=0.03)
    plt.legend(loc='best', ncol=3)
    
    # modify the legends
    new_labels = ['All','Photo + Post-GPP','Photo + Pheno','Photo','Pheno','Post-GPP']
    for ax in g.axes.flat:
        leg = ax.get_legend()
        if not leg is None: break
    # or legend may be on a figure
    if leg is None: leg = ax.legend
    
    for t, l in zip(leg.texts, new_labels): t.set_text(l)
    
    # entire figure title
    g.fig.suptitle('RMSE comparison ' + r"$\bf{WITH}$ " + 'different parameter sets' + ': ' + optim[6] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
                                
    # save figure
    g.fig.savefig(outdir + 'RMSE/RMSE_comparison_nee_optimization_different_param_set.png')
    plt.close()
    
    
    ## plotting NEE optimization MEAN RMSE comparison for various parameter sets
    ## ---------------------
    g = sns.catplot(x='Flux', y='RMSE', hue='optim', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','yellow','green','blue','orange','red'], height=6, 
                         legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.9,left=0.08,bottom=0.07)
    plt.legend(loc='best',ncol=2,title='Parameter sets')
    
    # modify the legends
    new_labels = ['All','Photo + Post-GPP','Photo + Pheno','Photo','Pheno','Post-GPP']
    for ax in g.axes.flat:
        leg = ax.get_legend()
        if not leg is None: break
    # or legend may be on a figure
    if leg is None: leg = ax.legend
    
    for t, l in zip(leg.texts, new_labels): t.set_text(l)
    
    # entire figure title
    g.fig.suptitle('Mean RMSE comparison ' + r"$\bf{WITH}$ " + 'different parameter sets' + '\n' + optim[6] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
    g.fig.text(y=0.02, x=0.45, verticalalignment='center', s='Flux', size=12, rotation=0)
                         
    # save figure
    g.fig.savefig(outdir + 'RMSE/Mean_RMSE_comparison_nee_optimization_different_param_set.png')
    plt.close()
    
    
    
    #%% plotting GPP+Reco optimization RMSE comparison for various parameter sets
    ## ---------------------
    ## ---------------------
if rmse_plot == True and [1,2,3,4,5,6,7] in df_rmse_all.optim.unique(): 
    df_rmse_all_plot = df_rmse_all.loc[df_rmse_all['optim'].isin([1,2,3,4,5,6,7])]
    df_rmse_all_plot = df_rmse_all_plot.loc[df_rmse_all_plot['Datatype'] == 'Ratio']
    
    df_rmse_all_plot = df_rmse_all_plot.drop('Datatype', 1)
    df_rmse_all_plot = df_rmse_all_plot.drop('optim_name', 1)
              
    g = sns.catplot(x='Flux', y='RMSE', hue='optim', col='Site', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','yellow','green','blue','orange','red'], height=6, 
                         col_wrap = col_wrap, col_order = site_names, legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.92,left=0.03)
    plt.legend(loc='best', ncol=3)
    
    # modify the legends
    new_labels = ['All','Photo + Post-GPP','Photo + Pheno','Photo','Pheno','Post-GPP']
    for ax in g.axes.flat:
        leg = ax.get_legend()
        if not leg is None: break
    # or legend may be on a figure
    if leg is None: leg = ax.legend
    
    for t, l in zip(leg.texts, new_labels): t.set_text(l)
    
    # entire figure title
    g.fig.suptitle('RMSE comparison ' + r"$\bf{WITH}$ " + 'different parameter sets' + ': ' + optim[0] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
                                
    # save figure
    g.fig.savefig(outdir + 'RMSE/RMSE_comparison_gppreco_optimization_different_param_set.png')
    plt.close()
    
    
    ## plotting GPP+Reco optimization MEAN RMSE comparison for various parameter sets
    ## ---------------------
    g = sns.catplot(x='Flux', y='RMSE', hue='optim', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','yellow','green','blue','orange','red'], height=6, 
                         legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(xlabel='', ylabel='')
    plt.subplots_adjust(top=0.9,left=0.08,bottom=0.07)
    plt.legend(loc='best',ncol=2,title='Parameter sets')
    
    # modify the legends
    new_labels = ['All','Photo + Post-GPP','Photo + Pheno','Photo','Pheno','Post-GPP']
    for ax in g.axes.flat:
        leg = ax.get_legend()
        if not leg is None: break
    # or legend may be on a figure
    if leg is None: leg = ax.legend
    
    for t, l in zip(leg.texts, new_labels): t.set_text(l)
    
    # entire figure title
    g.fig.suptitle('Mean RMSE comparison ' + r"$\bf{WITH}$ " + 'different parameter sets' + '\n' + optim[0] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
    g.fig.text(y=0.02, x=0.45, verticalalignment='center', s='Flux', size=12, rotation=0)
                                
    # save figure
    g.fig.savefig(outdir + 'RMSE/Mean_RMSE_comparison_gppreco_optimization_different_param_set.png')
    plt.close()
    
    
    #%% plotting GPP optimization RMSE comparison to test the contribution of Photosynthesis parameters
    ## ---------------------
    ## ---------------------
if rmse_plot == True and [15,16] in df_rmse_all.optim.unique(): 
    df_rmse_all_plot = df_rmse_all.loc[df_rmse_all['optim'].isin([15,16])]
    df_rmse_all_plot = df_rmse_all_plot.loc[df_rmse_all_plot['Datatype'] == 'Ratio']
    
    df_rmse_all_plot = df_rmse_all_plot.drop('Datatype', 1)
    df_rmse_all_plot = df_rmse_all_plot.drop('optim_name', 1)
              
    g = sns.catplot(x='optim', y='RMSE', hue='Flux', col='Site', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','blue','orange'], col_order = site_names, legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(ylabel='', xlabel='')
    plt.subplots_adjust(top=0.89,left=0.03,bottom=0.08)
    plt.legend(loc='best')
    
    # entire figure title
    g.fig.suptitle('RMSE comparison ' + r"$\bf{FOR}$ " + 'Photosynthesis parameters' + ': ' + optim[13] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
    g.fig.text(y=0.02, x=0.45, verticalalignment='center', s='Parameter set', size=12, rotation=0)
                                
    for ax in g.axes.flat:
        labels = ax.get_xticklabels() # get x labels
        labels[0] = 'Pheno + Photo' 
        labels[1] = 'Pheno' 
        ax.set_xticklabels(labels) # set new labels 
                            
    # save figure
    g.fig.savefig(outdir + 'RMSE/RMSE_comparison_gpp_optimization_photo_param.png')
    plt.close()
    
    ## plotting GPP optimization MEAN RMSE comparison to test the contribution of Photosynthesis parameters
    ## ---------------------
    g = sns.catplot(x='optim', y='RMSE', hue='Flux', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','blue','orange'], legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(ylabel='', xlabel='')
    plt.subplots_adjust(top=0.89,left=0.1,bottom=0.08)
    plt.legend(loc='best')
    
    # entire figure title
    g.fig.suptitle('Mean RMSE comparison ' + r"$\bf{FOR}$ " + 'Photosynthesis parameters' + '\n' + optim[13] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
    g.fig.text(y=0.02, x=0.45, verticalalignment='center', s='Parameter set', size=12, rotation=0)
                                
    for ax in g.axes.flat:
        labels = ax.get_xticklabels() # get x labels
        labels[0] = 'Pheno + Photo' 
        labels[1] = 'Pheno' 
        ax.set_xticklabels(labels) # set new labels 
                            
    # save figure
    g.fig.savefig(outdir + 'RMSE/Mean_RMSE_comparison_gpp_optimization_photo_param.png')
    plt.close()
    
    
    #%% plotting GPP optimization RMSE comparison to test the contribution of Phenology parameters
    ## ---------------------
    ## ---------------------
if rmse_plot == True and [15,17] in df_rmse_all.optim.unique(): 
    df_rmse_all_plot = df_rmse_all.loc[df_rmse_all['optim'].isin([15,17])]
    df_rmse_all_plot = df_rmse_all_plot.loc[df_rmse_all_plot['Datatype'] == 'Ratio']
    
    df_rmse_all_plot = df_rmse_all_plot.drop('Datatype', 1)
    df_rmse_all_plot = df_rmse_all_plot.drop('optim_name', 1)
              
    g = sns.catplot(x='optim', y='RMSE', hue='Flux', col='Site', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','blue','orange'], col_order = site_names, legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(ylabel='', xlabel='')
    plt.subplots_adjust(top=0.89,left=0.03,bottom=0.08)
    plt.legend(loc='best')
    
    # entire figure title
    g.fig.suptitle('RMSE comparison ' + r"$\bf{FOR}$ " + 'Phenology parameters' + ': ' + optim[13] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
    g.fig.text(y=0.02, x=0.45, verticalalignment='center', s='Parameter set', size=12, rotation=0)
                                
    for ax in g.axes.flat:
        labels = ax.get_xticklabels() # get x labels
        labels[0] = 'Pheno + Photo' 
        labels[1] = 'Photo' 
        ax.set_xticklabels(labels) # set new labels 
                            
    # save figure
    g.fig.savefig(outdir + 'RMSE/RMSE_comparison_gpp_optimization_phenology_param.png')
    plt.close()
    
    ## plotting GPP optimization MEAN RMSE comparison to test the contribution of Photosynthesis parameters
    ## ---------------------
    g = sns.catplot(x='optim', y='RMSE', hue='Flux', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','blue','orange'], legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(ylabel='', xlabel='')
    plt.subplots_adjust(top=0.89,left=0.1,bottom=0.08)
    plt.legend(loc='best')
    
    # entire figure title
    g.fig.suptitle('Mean RMSE comparison ' + r"$\bf{FOR}$ " + 'Phenology parameters' + '\n' + optim[13] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
    g.fig.text(y=0.02, x=0.45, verticalalignment='center', s='Parameter set', size=12, rotation=0)
                                
    for ax in g.axes.flat:
        labels = ax.get_xticklabels() # get x labels
        labels[0] = 'Pheno + Photo' 
        labels[1] = 'Photo' 
        ax.set_xticklabels(labels) # set new labels 
                            
    # save figure
    g.fig.savefig(outdir + 'RMSE/Mean_RMSE_comparison_gpp_optimization_phenology_param.png')
    plt.close()
    
    #%% plotting Reco optimization RMSE comparison to test the contribution of Post-GPP parameters
    ## ---------------------
    ## ---------------------
if rmse_plot == True and [18,19,20] in df_rmse_all.optim.unique(): 
    df_rmse_all_plot = df_rmse_all.loc[df_rmse_all['optim'].isin([18,19,20])]
    df_rmse_all_plot = df_rmse_all_plot.loc[df_rmse_all_plot['Datatype'] == 'Ratio']
    
    df_rmse_all_plot = df_rmse_all_plot.drop('Datatype', 1)
    df_rmse_all_plot = df_rmse_all_plot.drop('optim_name', 1)
              
    g = sns.catplot(x='optim', y='RMSE', hue='Flux', col='Site', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','blue','orange'], col_order = site_names, legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(ylabel='', xlabel='')
    plt.subplots_adjust(top=0.89,left=0.03,bottom=0.08)
    plt.legend(loc='best')
    
    # entire figure title
    g.fig.suptitle('RMSE comparison ' + r"$\bf{FOR}$ " + 'Post-GPP parameters' + ': ' + optim[16] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
    g.fig.text(y=0.02, x=0.45, verticalalignment='center', s='Parameter set', size=12, rotation=0)
                                
    for ax in g.axes.flat:
        labels = ax.get_xticklabels() # get x labels
        labels[0] = 'All' 
        labels[1] = 'Pheno + Photo'
        labels[2] = 'Post-GPP'
        ax.set_xticklabels(labels) # set new labels 
                            
    # save figure
    g.fig.savefig(outdir + 'RMSE/RMSE_comparison_reco_optimization_post-GPP_param.png')
    plt.close()

    ## plotting GPP optimization MEAN RMSE comparison to test the contribution of Photosynthesis parameters
    ## ---------------------
    g = sns.catplot(x='optim', y='RMSE', hue='Flux', data=df_rmse_all_plot, kind='bar',
                         palette=['grey','blue','orange'], legend=False)    
    g = g.map(plt.axhline, y=1, ls='--', c='red') # add a horizontal line at ratio = 1.0
    [plt.setp(ax.texts, text="") for ax in g.axes.flat]
    g.set_titles(col_template = '{col_name}')
    g.set(ylabel='', xlabel='')
    plt.subplots_adjust(top=0.89,left=0.1,bottom=0.08)
    plt.legend(loc='best')
    
    # entire figure title
    g.fig.suptitle('Mean RMSE comparison ' + r"$\bf{FOR}$ " + 'Post-GPP parameters' + '\n' + optim[16] + ' optimization')   
    # one ylabel
    g.fig.text(x=0.005, y=0.5, verticalalignment='center', s='RMSE posterior / RMSE prior', size=12, rotation=90)
    g.fig.text(y=0.02, x=0.45, verticalalignment='center', s='Parameter set', size=12, rotation=0)
                                
    for ax in g.axes.flat:
        labels = ax.get_xticklabels() # get x labels
        labels[0] = 'All' 
        labels[1] = 'Pheno + Photo'
        labels[2] = 'Post-GPP'
        ax.set_xticklabels(labels) # set new labels 
                            
    # save figure
    g.fig.savefig(outdir + 'RMSE/Mean_RMSE_comparison_reco_optimization_post-GPP_param.png')
    plt.close()



    #%% plotting manuscript figure #1
    ## ---------------------
    ## ---------------------     
if ms_plot == True:
    img1 = mpimg.imread(outdir + 'Scatter_plots/NEE_prior_P1_scatter_plots_for_NEE_optimization.png')
    img2 = mpimg.imread(outdir + 'Scatter_plots/NEE_prior_slope_inset_scatter_plots.png')
    img3 = mpimg.imread(outdir + 'Scatter_plots/NEE_P1_slope_inset_scatter_plots.png')
    img4 = mpimg.imread(outdir + 'Scatter_plots/NEE_prior_corr_inset_scatter_plots.png')
    img5 = mpimg.imread(outdir + 'Scatter_plots/NEE_P1_corr_inset_scatter_plots.png')
    
    fig, ax1 = plt.subplots()
    ax1 = plt.imshow(img1)
    plt.axis('off')
    # add the inset figures
    ax2 = fig.add_axes([0.33, 0.315, 0.15, 0.15])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = fig.add_axes([0.66, 0.315, 0.15, 0.15])
    ax3 = plt.imshow(img3)
    plt.axis('off')
    
    ax4 = fig.add_axes([0.17, 0.57, 0.15, 0.15])
    ax4 = plt.imshow(img4)
    plt.axis('off')
    ax5 = fig.add_axes([0.5, 0.57, 0.15, 0.15])
    ax5 = plt.imshow(img5)
    plt.axis('off')
    
#    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_1.png', dpi = 800, bbox_inches='tight')
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_1.png', dpi = 300, bbox_inches='tight')
    plt.close()
    
    
#    img1 = mpimg.imread(outdir + 'Scatter_plots/NEE_prior_P1_scatter_plots_for_NEE_optimization.png')
#    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_diff_sub_box_whisker.png')
#    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_anomaly_slope_sub_box_whisker.png')
#    img4 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_anomaly_corrcoef^2_sub_box_whisker.png')
#    
#    fig1 = plt.figure(constrained_layout=True)
#    gs1 = fig1.add_gridspec(nrows=2, ncols=3, left=0, right=1, wspace=0.005, hspace=0.1)
#    f9_ax1 = fig1.add_subplot(gs1[:-1, :])
#    f9_ax1 = plt.imshow(img1)
#    plt.axis('off')
#    
#    f9_ax2 = fig1.add_subplot(gs1[-1, 0])
#    f9_ax2 = plt.imshow(img2)
#    plt.axis('off')
#    f9_ax3 = fig1.add_subplot(gs1[-1, 1])
#    f9_ax3 = plt.imshow(img3)
#    plt.axis('off')
#    f9_ax4 = fig1.add_subplot(gs1[-1, -1])
#    f9_ax4 = plt.imshow(img4)
#    plt.axis('off')
#  
#    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_1_v2.png', dpi = 800, bbox_inches='tight')
#    plt.close()
    
    
    ## plotting manuscript figure #2
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Mean_monthly_fluxes/NEE_mean_monthly_fluxes_optim_NEE_P1.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_2.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #3
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_daily_P1_bias_variance_phase_sub_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_daily_P1_bias_variance_phase_sub_box_whisker.png')
    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_daily_P1_bias_variance_phase_sub_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0.005, hspace=0.005)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[0, 1])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = plt.subplot(gs[0, -1])
    ax3 = plt.imshow(img3)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_3.png', dpi = 1000, bbox_inches='tight')
    plt.close()
        
           
    ## plotting manuscript figure #4
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_only_daily_bias_variance_phase_sub_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_4.png', dpi = 1500, bbox_inches='tight')
    plt.close()
        
           
    ## plotting manuscript figure #S1
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_timeseries/NEE_optim_nee_iter3840_paramP1_Csink.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S1.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #S2
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_P1_anomaly_corrcoef_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_P1_anomaly_corrcoef_box_whisker.png')
    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_P1_anomaly_corrcoef_box_whisker.png')
    img4 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_P1_anomaly_slope_box_whisker.png')
    img5 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_P1_anomaly_slope_box_whisker.png')
    img6 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_P1_anomaly_slope_box_whisker.png')
    
    gs = gridspec.GridSpec(2, 3)
    gs.update(wspace=0.005, hspace=0.005)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[0, 1])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = plt.subplot(gs[0, 2])
    ax3 = plt.imshow(img3)
    plt.axis('off')
    ax4 = plt.subplot(gs[1, 0])
    ax4 = plt.imshow(img4)
    plt.axis('off')
    ax5 = plt.subplot(gs[1, 1])
    ax5 = plt.imshow(img5)
    plt.axis('off')
    ax6 = plt.subplot(gs[1, 2])
    ax6 = plt.imshow(img6)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S2.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #S3
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Mean_monthly_fluxes/Mean_monthly_fluxes_optim_NEE_P1.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S3.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #S4
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_P1_bias_variance_phase_sub_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_annual_P1_bias_variance_phase_sub_box_whisker.png')
    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_annual_P1_bias_variance_phase_sub_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0.005, hspace=0.005)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[0, 1])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = plt.subplot(gs[0, -1])
    ax3 = plt.imshow(img3)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S4.png', dpi = 1000, bbox_inches='tight')
    plt.close()
        
    
    ## plotting manuscript figure #S7
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Mean_monthly_fluxes/Mean_monthly_fluxes_optim_NEE_sitetype.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S7.png', dpi = 300, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #S6
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_anomaly_corrcoef_sub_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_annual_anomaly_corrcoef_sub_box_whisker.png')
    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_annual_anomaly_corrcoef_sub_box_whisker.png')
    img4 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_anomaly_slope_sub_box_whisker.png')
    img5 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_annual_anomaly_slope_sub_box_whisker.png')
    img6 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_annual_anomaly_slope_sub_box_whisker.png')
    
    gs = gridspec.GridSpec(2, 3)
    gs.update(wspace=0.005, hspace=0.005)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[0, 1])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = plt.subplot(gs[0, 2])
    ax3 = plt.imshow(img3)
    plt.axis('off')
    ax4 = plt.subplot(gs[1, 0])
    ax4 = plt.imshow(img4)
    plt.axis('off')
    ax5 = plt.subplot(gs[1, 1])
    ax5 = plt.imshow(img5)
    plt.axis('off')
    ax6 = plt.subplot(gs[1, 2])
    ax6 = plt.imshow(img6)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S6.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #S5
    ## ---------------------
    ## --------------------- 
    img1 = mpimg.imread(outdir + 'Scatter_plots/NEE_scatter_plots_for_NEE_optimization.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S5.png', dpi = 300, bbox_inches='tight')
    plt.close()
        
           
    ## plotting manuscript figure #S8
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_daily_bias_variance_phase_sub_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_daily_bias_variance_phase_sub_box_whisker.png')
    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_daily_bias_variance_phase_sub_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[0, 1])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = plt.subplot(gs[0, -1])
    ax3 = plt.imshow(img3)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S8.png', dpi = 1000, bbox_inches='tight')
    plt.close()
        
           
    ## plotting manuscript figure #S9
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Parameters/Parameters_deviation_with_uncertainty_reduction_NEE_opt.png')
    img2 = mpimg.imread(outdir + 'Parameters/Parameters_uncertainty_reduction_NEE_opt.png')
#    img1 = mpimg.imread(outdir + 'Parameters/Parameters_deviation_with_uncertainty_reduction_NEE_opt_subset.png')
#    img2 = mpimg.imread(outdir + 'Parameters/Parameters_uncertainty_reduction_NEE_opt_subset.png')
    
    gs = gridspec.GridSpec(2, 1)
    gs.update(wspace=0, hspace=0.005)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[1, 0])
    ax2 = plt.imshow(img2)
    plt.axis('off')
        
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S9.png', dpi = 600, bbox_inches='tight')
    plt.close()


    ## plotting manuscript figure #S10
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'bcorr/bcorr_Vcm_opt1_p1.png')
    img2 = mpimg.imread(outdir + 'bcorr/bcorr_Vcm_opt1_p2.png')
    img3 = mpimg.imread(outdir + 'bcorr/bcorr_Vcm_opt1_p3.png')
    img4 = mpimg.imread(outdir + 'bcorr/bcorr_Vcm_opt1_p4.png')
    img5 = mpimg.imread(outdir + 'bcorr/bcorr_Vcm_opt1_p5.png')
    img6 = mpimg.imread(outdir + 'bcorr/bcorr_Vcm_opt1_p6.png')
    img7 = mpimg.imread(outdir + 'bcorr/bcorr_Vcm_opt1_p7.png')
    
    fig9 = plt.figure(constrained_layout=False)
    gs1 = fig9.add_gridspec(nrows=2, ncols=1, left=0, right=0.35, wspace=0.005, hspace=0.1)
    f9_ax1 = fig9.add_subplot(gs1[0, -1])
    f9_ax1.set_title('P1',fontsize=6)
    f9_ax1 = plt.imshow(img1)
    plt.axis('off')
    f9_ax2 = fig9.add_subplot(gs1[-1, -1])
    f9_ax2.set_title('P2',fontsize=6)
    f9_ax2 = plt.imshow(img2)
    plt.axis('off')
    
    gs2 = fig9.add_gridspec(nrows=2, ncols=1, left=0.35, right=0.7, wspace=0.005, hspace=0.1)
    f9_ax3 = fig9.add_subplot(gs2[0, -1])
    f9_ax3.set_title('P3',fontsize=6)
    f9_ax3 = plt.imshow(img3)
    plt.axis('off')
    f9_ax4 = fig9.add_subplot(gs2[-1, -1])
    f9_ax4.set_title('P4',fontsize=6)
    f9_ax4 = plt.imshow(img4)
    plt.axis('off')
    
    gs3 = fig9.add_gridspec(nrows=3, ncols=1, left=0.7, right=1, wspace=0.005, hspace=0.1)
    f9_ax5 = fig9.add_subplot(gs3[0, :])
    f9_ax5.set_title('P5',fontsize=6)
    f9_ax5 = plt.imshow(img5)
    plt.axis('off')
    f9_ax6 = fig9.add_subplot(gs3[1, :])
    f9_ax6.set_title('P6',fontsize=6)
    f9_ax6 = plt.imshow(img6)
    plt.axis('off')
    f9_ax7 = fig9.add_subplot(gs3[-1, :])
    f9_ax7.set_title('P7',fontsize=6)
    f9_ax7 = plt.imshow(img7)
    plt.axis('off')
    
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S10.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #S11
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Parameters/Parameters_all_optim_nee_iter3840_paramP1_Csink.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_S11.png', dpi = 1500, bbox_inches='tight')
    plt.close()
    
    
        
    
    
    
#%% plotting extra figures    
    ## plotting manuscript figure #E1
    ## ---------------------
    ## ---------------------  
    img1 = mpimg.imread(outdir + 'Scatter_plots/NEE_prior_P1_scatter_plots_for_NEE_optimization.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
  
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E1.png', dpi = 800, bbox_inches='tight')
    plt.close()
        
    
    
    ## plotting manuscript figure #E2
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_monthly_bias_variance_phase_sub_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_monthly_bias_variance_phase_sub_box_whisker.png')
    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_monthly_bias_variance_phase_sub_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0.005, hspace=0.005)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[0, 1])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = plt.subplot(gs[0, -1])
    ax3 = plt.imshow(img3)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E2.png', dpi = 1000, bbox_inches='tight')
    plt.close()
        
           
    ## plotting manuscript figure #E3
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_bias_variance_phase_sub_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_annual_bias_variance_phase_sub_box_whisker.png')
    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_annual_bias_variance_phase_sub_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0.005, hspace=0.005)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[0, 1])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = plt.subplot(gs[0, -1])
    ax3 = plt.imshow(img3)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E3.png', dpi = 1000, bbox_inches='tight')
    plt.close()
        
           
    ## plotting manuscript figure #E4
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Mean_monthly_fluxes/NEE_mean_monthly_fluxes_optim_NEE_all.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E4.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #E5
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_diff_sub_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_anomaly_slope_sub_box_whisker.png')
    
    gs = gridspec.GridSpec(1,2)
    gs.update(wspace=0.005, hspace=0.05)
    ax1 = plt.subplot(gs[0,0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[0,1])
    ax2 = plt.imshow(img2)
    plt.axis('off')
  
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E5.png', dpi = 800, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #E6
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_anomaly_corrcoef_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_annual_anomaly_corrcoef_box_whisker.png')
    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_annual_anomaly_corrcoef_box_whisker.png')
    
    gs = gridspec.GridSpec(3, 1)
    gs.update(wspace=0.005, hspace=0.005)
    ax1 = plt.subplot(gs[0,0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[1,0])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = plt.subplot(gs[-1,0])
    ax3 = plt.imshow(img3)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E6.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #E7
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_anomaly_bias_variance_phase_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E7.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #E8
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_annual_anomaly_bias_variance_phase_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E8.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #E9
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_annual_anomaly_bias_variance_phase_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E9.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #E10
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Mean_monthly_fluxes/Mean_monthly_fluxes_optim_GPP+Reco_all.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E10.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #E11
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Parameters/Parameters_deviation_with_uncertainty_reduction.png')
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.005, hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E11.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    ## plotting manuscript figure #E12
    ## ---------------------
    ## ---------------------     
    img1 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/NEE_flux_annual_anomaly_bias_variance_phase_sub_fraction_box_whisker.png')
    img2 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/GPP_flux_annual_anomaly_bias_variance_phase_sub_fraction_box_whisker.png')
    img3 = mpimg.imread(outdir + 'Flux_annual_anomaly/box_whisker_plots/Reco_flux_annual_anomaly_bias_variance_phase_sub_fraction_box_whisker.png')
    
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0.005, hspace=0.005)
    ax1 = plt.subplot(gs[0, 0])
    ax1 = plt.imshow(img1)
    plt.axis('off')
    ax2 = plt.subplot(gs[0, 1])
    ax2 = plt.imshow(img2)
    plt.axis('off')
    ax3 = plt.subplot(gs[0, -1])
    ax3 = plt.imshow(img3)
    plt.axis('off')
            
    plt.savefig(outdir + 'Manuscript_figures/' + 'figure_E12.png', dpi = 1000, bbox_inches='tight')
    plt.close()
    
    
    
#%% plotting manuscript figure #S
    ## ---------------------
    ## ---------------------    



    

