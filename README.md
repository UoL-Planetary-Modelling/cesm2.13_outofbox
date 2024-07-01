# WACCMX_metal_ion_transport_SMax
- Source code and required input files for the WACCMX Solar Maximum run with metal_ion_transport (Case name: SMax_3M_FX2000_f19f19mg16)

# Contents
- SourceMods/src.cam/ folder
- Chemistry input files
- User_nl_cam
- env_build.xml
- env_run.xml
- SMax.run submission script

# Case Directory located at:
"/nobackup/sestay/cesm2/cases/SMax_3M_FX2000_f19f19mg16/"

# File archive located at:
/nobackup/sestay/cesm2/archive/SMax_3M_FX2000_f19f19mg16/

# User_nl_cam
- SourceMods, and env_build/run are the same for SMin and SMax runs, differences in user_nl_cam:

`&solar_data_opts
 solar_data_type                = 'FIXED'
solar_data_ymd         = 20000101
 solar_htng_spctrl_scl          = .true.
 solar_irrad_data_file          = '/home/home02/sestay/WACCM_Input_Files/Smax_spectral_irradiance_cycle22.nc'
solar_parms_data_file          = '/home/home02/sestay/WACCM_Input_Files/Smax_params.nc'`
 
# Solar Input Files
- Solar input files linked in user_nl_cam are located at: /home/home02/sestay/WACCM_Input_Files/
- Files were too big to upload, but code to reproduce the files is in SpE_Identification_Notebooks repo 
    - Link here: https://github.com/tashaay/SpE_Identification_Notebooks/blob/main/Solar%20Irradiance.ipynb
