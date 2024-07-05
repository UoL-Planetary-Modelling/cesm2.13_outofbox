**WACCM-X simulation with electrodynamical ion transport for Fe, Mg and Na under solar maximum conditions**

**Configuration**
If you haven't already, follow the steps in configuration from the FX2000 file in the main branch


**Set up new case and change to the case directory:**
`/home/home01/earfw/release_cesm2_1_3/cime/scripts/create_newcase --case /nobackup/$USER/cesm2/cases/f.e21.FWmadSD.f09_f09_mg17.2005-2020.002 --compset FX2000 --res f19_f19_mg16 --machine arc4 --run-unsupported`
`cd /nobackup/$USER/cesm2/cases/f.e21.FWmadSD.f09_f09_mg17.2005-2020.002`


**Change number of Cores**
`./xmlchange NTASKS=80`


**Case Setup**
`./case.setup`

**Case modifications**
Copy in user_nl_cam, sourcemods, env_build and env_run files from github to case directory (clone the online repo to the remote server (arc4) and then copy the files into the case directory from there). 

`git clone https://github.com/UoL-Planetary-Modelling/cesm2.13/ [folder path to store repo on remote server]`
N.B. this clones the entire repo

`git checkout f.e21.FWmadSD.f09_f09_mg17.2005-2020.002` - swap to the right branch

`cp [folder path to store repo on remote server]/SourceMods/src.cam/* /nobackup/$USER/cesm2/cases/f.e21.FWmadSD.f09_f09_mg17.2005-2020.002`
`cp [folder path to store repo on remote server]/env_run.xml /nobackup/$USER/cesm2/cases/f.e21.FWmadSD.f09_f09_mg17.2005-2020.002`
`cp [folder path to store repo on remote server]/env_build.xml /nobackup/$USER/cesm2/cases/f.e21.FWmadSD.f09_f09_mg17.2005-2020.002`
`cp [folder path to store repo on remote server]/SMin.run /nobackup/$USER/cesm2/cases/f.e21.FWmadSD.f09_f09_mg17.2005-2020.002`


Ensure links to input files in user_nl_cam are pointing to the right files. For the run I did, the files were stored here: /nobackup/sestay/WACCM_Input_Files/ but I have copied them to here /resstore/b0243/Data/Input_Files_f.e21.FWmadSD.f09_f09_mg17.2005-2020.002/ to avoid deletion by arc /nobackup/

  solar_irrad_data_file = '/nobackup/sestay/WACCM_Input_Files/SolarForcingNRLSSI2_daily_s18820101_e20201231_extended_c20210414.nc'
  solar_parms_data_file = '/nobackup/sestay/WACCM_Input_Files/SolarGeomag_s19491230_e20210411_c20210414.nc'
  epp_all_filepath      = '/nobackup/sestay/WACCM_Input_Files/SolarForcingCCMI2022REFD1_19491230-20200101_sumEPP_extended_c20210414.nc'


**Build the model**
`./case.build`
`./case.build --skip-provenance-check`

**Submit the run**
Submit the job: qsub -cwd WACCM_D.run
