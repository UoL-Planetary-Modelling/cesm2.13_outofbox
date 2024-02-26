# Steps to run FWmadSD in CESM2_1_3 on ARC4
## Wuhu Feng, NCAS, University of Leeds 14 October 2023

1) add the following lines in your ~/.bashrc, then logout arc4 then relogin (only do this once)

>
    if [ "$SERVICE_NAME" = "arc4" ]; then
    export CIME_MODEL=cesm
    export CCSMUSER=$USER
    if [ ! -d /nobackup/$USER/cesm2_inputdata ]; then
    ln -s /nobackup/earfw/cesm2_inputdata/nobackup/$USER/cesm2_inputdata
    fi
    # Load Modules
    if [ -r /nobackup/cemac/cemac.sh ] ; then
        . /nobackup/cemac/cemac.sh
    fi
    #intel
    module load netcdf
    module load python3
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



2) create a new case (for example specified-version of WACCM-D at 1 degree).

>
    /home/home01/earfw/release_cesm2_1_3/cime/scripts/create_newcase  --case /nobackup/$USER/cesm2/cases/f.e213.FWmadSD.f09_f09_mg17 --compset FWmadSD --res f09_f09_mg17 --machine arc4 --run-unsupported 

3) access the case directory:
>
    cd /nobackup/$USER/cesm2/cases/f.e213.FWmadSD.f09_f09_mg17

4) Note the default simulation uses 160 cores, this can be changed either by xmlchange command or manually changes for the env_mach_pes.xml, for example, 1 node (40 cores), 2 nodes (80 cores).
Here just gives one example using 1 dedicate planet node:

Method 1: ./xmlchange NTASKS=-1

Method 2: vi env_mach_pes.xml, then change -4 to -1.

5) set up the case:
>
    ./case.setup

6) Adding the follow lines in the user_nl_cam to make sure the MERRA2 reanalysis in the correct directory (only 2005 completed, but I will add other reanalysis years once the downloads are completed).
>
    &metdata_nl
     met_data_file          = '2005/MERRA2_0.9x1.25_20050101.nc'
     met_data_path          = '/nobackup/earfw/cesm2_inputdata//atm/cam/met/MERRA2/0.9x1.25'
     met_filenames_list     = '/nobackup/earfw/cesm2_inputdata//atm/cam/met/MERRA2/0.9x1.25/merra2_2005afterwards.txt'
    /

7) build the case:
>
    ./case.build

8) check if any input files are missing, note that you can not download the missing input files to Wuhu's directory (/nobackup/earfw/cesm2_inputdata):
>
     ./check_input_data --download

9) create the following jobscript for example WACCMD.run (to be submitted anywhere), note the case directory in the jobscript:

>
    #!/bin/csh -f
    #===============================================================================
    #Wuhu cesm2_2_0_wuhu
    #===============================================================================
    ################################################################################
    #Leeds University arc4: WUHU FENG
    #################################################################################
    #$ -V
    #$ -cwd
    #$ -l h_rt=48:00:00
    #$ -j y
    #$ -P planet
    #$ -l h_vmem=4G
    #$ -l np=40

    cd /nobackup/$USER/cesm2/cases/f.e213.FWmadSD.f09_f09_mg17
    python case.submit

10) submit the job :
>
    qsub -cwd WACCMD.run


11) If 5-day runs successful, then do the continous run following NCAR CESM instruction....
