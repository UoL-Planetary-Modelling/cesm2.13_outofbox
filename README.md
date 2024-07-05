# WACCMX_metal_ion_transport_SMin
- Source code and required input files for the WACCMX Solar Minimum run with metal_ion_transport


# Contents
- HOW_TO_RUN file
- SourceMods/src.cam/ folder
- Chemistry input files
- User_nl_cam
- env_build.xml
- env_run.xml
- SMin.run submission script

# Case Directory located at:
"/nobackup/sestay/cesm2/cases/SMin_3M_FX2000_f19f19mg16/"

# File archive located at:
/nobackup/sestay/cesm2/archive/SMin_3M_FX2000_f19f19mg16/

# Solar min conditions are:
- f10.7=70, f10.7a=70, kp=0.3, ap=1.8

# User_nl_cam
- SourceMods, and env_build/run are the same for SMin and SMax runs, differences in user_nl_cam:

```&solar_data_opts

 solar_data_type                = 'FIXED'
 
solar_data_ymd         = 20000101

 solar_htng_spctrl_scl          = .true.
 
 solar_irrad_data_file          = '/home/home02/sestay/WACCM_Input_Files/Smin_spectral_irradiance_cycle21.nc'
 
solar_parms_data_file          = '/home/home02/sestay/WACCM_Input_Files/Smin_params.nc'```
 
# Solar Input Files
- Solar input files linked in user_nl_cam are located at: /home/home02/sestay/WACCM_Input_Files/
- Files were too big to upload, but code to reproduce the files is in SpE_Identification_Notebooks repo 
    - Link here: https://github.com/tashaay/SpE_Identification_Notebooks/blob/main/Solar%20Irradiance.ipynb

# .bashrc file
- Add the lines to your .bashrc file to configure the ESMF libraries and modules used