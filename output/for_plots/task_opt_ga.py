import numpy as np
import os

path_flux = "/home/users/vbastri/FLUX/sites15.txt"
dtype = [("ID", "S6"), ("pft", "i2"), ("lat", "f8"), ("lon", "f8"), ("years", "S9"), ("veg", "S64"), ("soilt", "S65")]
sites = np.loadtxt(path_flux, delimiter = "\t", dtype=dtype)

start = """<?xml version="1.0"?>

<xml>
    <task id = "optimize"
          norc = "8"
          orchidee = "TRUNK"
          minimizer = "GA"
          maxiter = "25"
          path_output = "PFT%(pft)s"
          param_opt = "A1, B1, CHOISNEL_RSOL_CSTE, DEPTH_MAX_H, FRAC_GROWTHRESP, HCRIT_LITTER, HYDROL_HUMCSTE, KSOILC, LAI_MAX, LAI_MAX_TO_HAPPY, LEAFAGECRIT, MAINT_RESP_COEFF, MAINT_RESP_SLOPE_C, MOIST_COEFF, MOISTCONT_MIN, RESIDENCE_TIME, SENESCENCE_TEMP_C, SLA, SOIL_Q10, TAU_LEAFINIT, VCMAX25, Z_DECOMP, Z0_OVER_HEIGHT" />

    <minimizer id = "GA" population = "32" />
"""

flux = """
    <site id = "%(id)s"
          pft = "%(pft)s"
          type = "flux"
          obs_path = "/home/users/vbastri/FLUX/%(pft)s/%(id)s/%(id)s_flux_%(years)s.nc"
          obs_name = "NEEt, Qle"
          obs_time = "time_counter"
          obs_filter = "flatten, mean, smooth"
          obs_coef = "60*60*24, 1"
          sim_name = "NEEt, fluxlat"
          sim_filter = "flatten, smooth"
          sim_coef = "60*60*24, 1"
          forcing_path = "/home/users/vbastri/FLUX/%(pft)s/%(id)s/%(id)s%(years)s.nc"
          forcing_time = "time"
          sechiba_path = "/home/satellites7/vbastri/OPT/OPT4661/SPINUP/%(id)s/%(id)s_sechiba.nc"
          stomate_path = "/home/satellites7/vbastri/OPT/OPT4661/SPINUP/%(id)s/%(id)s_stomate.nc"
          sechiba_step = "86400" />
"""

end = """
</xml>
"""

for pft in set(sorted(sites["pft"])):
    task = start % dict(pft = pft)
    for site in sites[sites["pft"] == pft]:
        meta = dict(id = site["ID"], pft = site["pft"], years = site["years"], veg = site["veg"], soilt = site["soilt"])
        task += flux % meta
    task += end
    with open("pft%s.xml" % pft, "w") as stream: stream.write(task)

