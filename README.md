# cesm2.13_outofbox
How to run a 5 day 'out of the box' test simulation using CESM2.13 on arc4 - model version ported by Wuhu Feng (w.feng@ncas.ac.uk)

## Configuration 
- First you need to configure the environment for running the model - i.e. set up environment variables, load modules, link input data directories. You can do this by adding this information to your .bashrc file, which provides settings to be executed when a user starts an interactive shell session
`vi ~/.bashrc`
- Press 'i' to enter insert mode and paste in the following:
```
################ WUHU ADDED SECTION FOR USING WACCM ##################
if [ "$SERVICE_NAME" = "arc4" ]; then
export CIME_MODEL=cesm
export CCSMUSER=$USER
if [ ! -d /nobackup/$USER/cesm2_inputdata ]; then
ln -s /nobackup/earfw/cesm2_inputdata /nobackup/$USER/cesm2_inputdata
fi
# Load Modules
if [ -r /nobackup/cemac/cemac.sh ] ; then
   . /nobackup/cemac/cemac.sh
fi
#intel
module load netcdf
module load python/2.7.16
module load cmake
module load mkl
module load hdf5
export ESMF_DIR=/resstore/b0154/Data/earfw/CAM_6_3_110/esmf_8_4_2_arc4
export ESMFMKFILE=/resstore/b0154/Data/earfw/CAM_6_3_110/esmf_8_4_2_arc4/lib/libO/Linux.intel.64.openmpi.default/esmf.mk
export ESMF_COMM=openmpi
export ESMF_COMPILER=intel
export ESMF_ABI=64
export ESMF_OS=Linux
export NETCDF=$NETCDF_HOME
export NETCDF_PATH=$NETCDF_HOME
export DIN_LOC_ROOT=/nobackup/earfw/cesm2_inputdata
fi
#####################################################################
```
- Press ':wq' to save and exit
- Log out and back in again to source the .bashrc
- You can check what modules you currently have running using:
`module list`

## Python Environment
- If you have set up any specific python environments, deactivate these, as these can cause issues:
`conda deactivate`

## Create Case
- Create a new folder to host your case directories:
`mkdir /nobackup/$USER/cesm2/cases`

-	Set up a new case – this is an out of the box 5 day test run using cesm 2.13 compset FX2000 (WACCM-X perpetual year 2000) with 1.9x2.5degree resultion. You can change the case name (CESM213_FX2000_f19_f19_mg16_arc4) to something descriptive
`/home/home01/earfw/release_cesm2_1_3/cime/scripts/create_newcase --case /nobackup/$USER/cesm2/cases/CESM213_FX2000_f19_f19_mg16_arc4 --compset FX2000 --res f19_f19_mg16 --machine arc4 --run-unsupported`

- Change to the case directory:
`cd /nobackup/$USER/cesm2/cases/CESM213_FX2000_f19_f19_mg16_arc4`

## Case Setup
`./case.setup`

## Build the model
./case.build --skip-provenance-check

## Submit the run
- Create a submission script in the case directory
`vi /nobackup/$USER/cesm2/cases/CESM213_FX2000_f19_f19_mg16_arc4/submission_script.sh`
- Press 'i' to enter insert mode and paste in the following:
```
######################################################################
#!/bin/csh -f
#Run with current environment (-V) and in the current directory (-cwd)
#$ -V
#$ -cwd
#Request some time- min 15 mins - max 48 hours
#$ -l h_rt=08:00:00
#Combine the standard error and standard output into one file
#$ -j y
#Sets project if needed - change to relevant quue or remove this line to use the standard arc queue
#$ -P feps-cpu

#Specifies a job for parallel programs using MPI. Assigns whole compute nodes. in np=x, x is the number of processes
# Change number of cores if relevant. Standard is two (np=80) 
#$ -l np=40
#Change to case directory location
cd /home/home02/sestay/cesm_prep/test_FX2000_f19f19 
python case.submit
######################################################################
```
- Press ':wq' to save and exit
- N.B. you can change the line `#$ -l np=40` to `#$ -pe ib 40`. The latter will use the available shared resources from arc4 which can reduce the queueing time, but sometimes the model occasional crashing due to the shared available memory.

- Submit the job:
`qsub -cwd submission_script.sh` 



