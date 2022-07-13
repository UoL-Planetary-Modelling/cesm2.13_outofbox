# WACCMX_metal_ion_transport_SMax
- Source code and required input files for the WACCMX Solar Maximum run with metal_ion_transport
- SourceMods, input files and env_build are the same for SMin and SMax runs, differences in user_nl_cam 
      - SMin date for solar data files = 2008-11-01
              - Solar parameters are:
                      - f107 = 65.6
                      - f107a = 67.487656
                      - kp = 0.625
                      - ap = 2.75
      - SMax date for solar data files = 2001-11-01
                      - f107 = 232
                      - f107a = 215.9605
                      - kp = 3.3875
                      - ap = 24.125
                   
## Contents
- SourceMods/src.cam/ folder
- Chemistry input files
- User_nl_cam
- env_build.xml

## Case Directory located at:
/home/home02/sestay/cesm_prep/cesm213_MIon_SMax_FX2000_f19f19mg16_2/

## bld/ and run/ folders located at:
/nobackup/sestay/cesm_sims/cesm213_MIon_SMax_FX2000_f19f19mg16_2/

## File archive located at:
/nobackup/sestay/cesm_sims/archive/cesm213_MIon_SMax_FX2000_f19f19mg16_2/

## Input files located at:
/home/home02/sestay/WACCM_Input_Files/
