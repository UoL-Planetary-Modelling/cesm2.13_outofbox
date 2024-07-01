
if [ "$SERVICE_NAME" = "arc4" ]; then
export CIME_MODEL=cesm
export CCSMUSER=$USER
if [ ! -d /nobackup/$USER/cesm2_inputdata ]; then
ln -s /nobackup/earfw/cesm2_inputdata /nobackup/$USER/cesm2_inputdata
fi
#intel
module purge
module load licenses
module load sge
module load user
module load intel
module load openmpi
module load netcdf/4.6.3
module load cmake/3.15.1
module load python/2.7.16
module load mkl
module load hdf5

export ESMF_DIR=/nobackup/earfw/CESM2/esmf_arc4_03102022
export ESMFMKFILE=/nobackup/earfw/CESM2/esmf_arc4_03102022/lib/libO/Linux.intel.64.openmpi.default/esmf.mk
export ESMF_COMM=openmpi
export ESMF_COMPILER=intel
export ESMF_ABI=64
export ESMF_OS=Linux
export LD_LIBRARY_PATH=/nobackup/earfw/CESM2/esmf_arc4_03102022/lib:/nobackup/earfw/CESM2/esmf_arc4_03102022/lib/libO/Linux.intel.64.openmpi.default:${LD_LIBRARY_PATH}
export NETCDF=$NETCDF_HOME
export NETCDF_PATH=$NETCDF_HOME
export DIN_LOC_ROOT=/nobackup/earfw/cesm2_inputdata
fi
