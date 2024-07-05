**WACCM-X simulation with electrodynamical ion transport for Fe, Mg and Na under solar maximum conditions**

**Configuration**

If you haven't already, follow the steps in configuration from the FX2000 file in the main branch


**Create Case**

Set up new case:
`/home/home01/earfw/release_cesm2_1_3/cime/scripts/create_newcase --case /nobackup/$USER/cesm2/cases/SMax_3M_FX2000_f19f19mg16 --compset FX2000 --res f19_f19_mg16 --machine arc4 --run-unsupported`

Change to the case directory: 
`cd /nobackup/$USER/cesm2/cases/SMax_3M_FX2000_f19f19mg16`


**Change number of Cores**

`./xmlchange NTASKS=80`


**Case Setup**

`./case.setup`


**Case modifications**

Next you want to copy in user_nl_cam, sourcemods, env_build and env_run and SMax.run files from github to your case directory. 
You can clone the online repo to the remote server (arc4) and then copy the files into the case directory from there. 

`git clone https://github.com/UoL-Planetary-Modelling/cesm2.13/ [folder path to store repo on remote server]`

N.B. this clones the entire repo

`git checkout FX2000_SMax` - swap to the right branch

`cp [folder path to store repo on remote server]/SourceMods/src.cam/* /nobackup/$USER/cesm2/cases/SMax_3M_FX2000_f19f19mg16`

`cp [folder path to store repo on remote server]/env_run.xml /nobackup/$USER/cesm2/cases/SMax_3M_FX2000_f19f19mg16`

`cp [folder path to store repo on remote server]/env_build.xml /nobackup/$USER/cesm2/cases/SMax_3M_FX2000_f19f19mg16`

`cp [folder path to store repo on remote server]/SMin.run /nobackup/$USER/cesm2/cases/SMax_3M_FX2000_f19f19mg16`


**Build the model**

`./case.build`
`./case.build --skip-provenance-check`


**Submit the run**
Submit the job: qsub -cwd SMax.run
