#!/bin/bash -x
#SBATCH -J labseaHRspinup 
#SBATCH -o labseaHRspinup.%j.out
#SBATCH -e labseaHRspinup.%j.err
#SBATCH -t 96:00:00
#SBATCH -N 7 
#SBATCH -n 180

#--- 0.load modules ------
#module purge
##module load intel openmpi netcdf-fortran
#module load intel/2023.1.0 openmpi4/4.1.5 phdf5/1.14.1 netcdf-fortran/4.6.0 netcdf/4.9.0 prun
#echo $LD_LIBRARY_PATH

#ulimit -s hard
#ulimit -u hard


#---- set variables ------
# note: for nprocs, take ntiles - length(blanklist)
nprocs_hr=180

jobfile=submit_HR_spinup.sh

#--- set dir ------------
rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
datadir=/home/shoshi/MITgcm_obsfit/lab_sea12/
scratchdir=/scratch/shoshi/labsea_MG_12/

builddir=${rootdir}/build_hi_fwd


workdir_hi=${scratchdir}/run_hi_July24

if [ ! -d $workdir_hi ]; then
  mkdir -p $workdir_hi;
fi

cd $workdir_hi
mkdir diags

#--- 6. NAMELISTS ---------
#ln -s ${datadir}/input_cal/* .
cp ${rootdir}/input_adhi/* .
ln -s ${datadir}/input_binaries_hires/bathy_cleaned.bin .
ln -s ${datadir}/input_binaries_hires/LevCli_temp_480x384_linearv5.labsea1979 .
ln -s ${datadir}/input_binaries_hires/LevCli_salt_480x384_linearv5.labsea1979 .
ln -s ${datadir}/input_binaries_exf/* .
ln -s ${datadir}/input_binaries_hires/ones_64b.bin .
ln -s ${datadir}/input_binaries_hires/ARGO_WO_2024_PFL_D_labsea_splitcost.nc .
ln -s ${datadir}/input_weights_hires/*_jra3q_weights_Jan2024_64b_SMOOTHED_removeboundary.bin .
ln -s ${datadir}/input_weights_hires/*fromASTE_*.bin .
ln -s ${datadir}/input_binaries_hires/diffkr_r4_HR.bin .
#cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup_seaice.0000103680.meta . 
#cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup_seaice.0000103680.data .
#cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup.0000103680.meta . 
#cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup.0000103680.data . 
cp /scratch/shoshi/labsea_MG_12//run_hi_6mo_Jan24_to_June24/pickup_seaice.0000155520.meta . 
cp /scratch/shoshi/labsea_MG_12//run_hi_6mo_Jan24_to_June24/pickup_seaice.0000155520.data .
cp /scratch/shoshi/labsea_MG_12//run_hi_6mo_Jan24_to_June24/pickup.0000155520.meta . 
cp /scratch/shoshi/labsea_MG_12//run_hi_6mo_Jan24_to_June24/pickup.0000155520.data . 
cp ${datadir}/input_binaries_hires/smooth2Dscales001_4x.bin ./smooth2Dscales000.data
cp ${datadir}/input_binaries_hires/smooth2Dscales001_4x.bin ./smooth2Dscales001.data
cp ${datadir}/input_binaries_hires/smooth3DscalesH001_4x.bin ./smooth3DscalesH001.data
cp ${datadir}/input_binaries_hires/smooth3DscalesZ001_4x.bin ./smooth3DscalesZ001
cp ${datadir}/input_binaries_hires/smooth*operator* .
cp ${datadir}/input_binaries_hires/smooth2Doperator001.meta ./smooth2Doperator000.meta
cp ${datadir}/input_binaries_hires/smooth2Doperator001.data ./smooth2Doperator000.data
cp ${datadir}/input_binaries_hires/smooth2Dnorm001.meta ./smooth2Dnorm000.meta
cp ${datadir}/input_binaries_hires/smooth2Dnorm001.data ./smooth2Dnorm000.data
cp ${datadir}/input_binaries_hires/smooth*norm* .
mkdir jra3q
ln -s /scratch/shared/jra3q/jra3q_*_2024 ./jra3q/

#-- swap out data.ctrl and copy high-res adjustments
sed -i -e 's/'"doinitxx = .FALSE."'/'"doinitxx = .TRUE."'/g' data.ctrl
sed -i -e 's/'"doInitXX = .FALSE."'/'"doInitXX = .TRUE."'/g' data.ctrl
sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
sed -i -e 's/'"doMainUnpack = .TRUE."'/'"doMainUnpack = .FALSE."'/g' data.ctrl

rm data
mv data_spinup data
rm data.pkg
mv data.pkg_spinup data.pkg
rm data.cal
mv data.cal_spinup data.cal
rm data.diagnostics
cp /scratch/shoshi/labsea_MG_12/run_hi_1yr_jra3q_B/data.diagnostics .

\rm data.optim
cat > data.optim <<EOF
 &OPTIM
 optimcycle=0,
 /
EOF


#--- 7. executable --------
cp -p $builddir/mitgcmuv .


set -x
date > run.MITGCM.timing
mpiexec -np ${nprocs_hr} ./mitgcmuv > stdout
date >> run.MITGCM.timing
cd ..

