#!/bin/bash -x
#SBATCH -J labseaHRswot_adj
#SBATCH -o labseaHRswot_adj.%j.out
#SBATCH -e labseaHRswot_adj.%j.err
#SBATCH -t 6:00:00
#SBATCH -N 6 
#SBATCH -n 103
#SBATCH --mail-user=sreich@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


iter=0
optimext=0000
nprocs_hr=103

#--- 2.set dir ------------
rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
datadir=/home/shoshi/MITgcm_obsfit/lab_sea12/
scratchdir=/scratch/shoshi/labsea_MG_12/assim_swot_HR/
#dirrun_pup=${scratchdir}/run_adlo_packunpack/

hr_dir_iter0=${scratchdir}/run_adhi_it0000

builddir=${rootdir}/build_adhi_2lev_seaice_update_mpi_32x32

#echo $(printf "%04d" $((iter)))  
#source $(conda info --base)/etc/profile.d/conda.sh
#conda activate py38
#python3 ${rootdir}/mappings/interp_xx_lores_to_hires_itX_v2.py $(printf "%04d" $((iter))) $dirrun_pup
#conda deactivate
#

ext2=$(printf "%04d" $iter)


workdir_hi=${scratchdir}/run_adhi_it${ext2}

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
ln -s ${datadir}/input_binaries_hires/swot_obsfit_cycles_9thru11_labsea_L3v3.nc .
ln -s ${datadir}/input_binaries_hires/swot_obsfit_cycles_9thru11_labsea_L3v3_PROFILES.nc .
ln -s ${datadir}/input_binaries_hires/rads_20240101_20240308.nc .
ln -s ${datadir}/input_binaries_hires/rads_* .
ln -s ${datadir}/input_binaries_lores/slaerr_03m.bin .
ln -s ${datadir}/input_weights_hires/*_jra3q_weights_Jan2024_64b_SMOOTHED_removeboundary.bin .
ln -s ${datadir}/input_weights_hires/*fromASTE_*.bin .
ln -s ${datadir}/input_binaries_hires/diffkr_r4_HR.bin .
cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup_seaice.0000103680.meta pickup_seaice.0000000001.meta
cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup_seaice.0000103680.data pickup_seaice.0000000001.data
cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup.0000103680.meta pickup.0000000001.meta
cp /scratch/shoshi/labsea_MG_12//run_hi_1yr_jra3q_B/pickup.0000103680.data pickup.0000000001.data
#cp ${datadir}/input_binaries_hires/pickup_seaice.0000103680.meta pickup_seaice.0000000001.meta
#cp ${datadir}/input_binaries_hires/pickup_seaice.0000103680.data pickup_seaice.0000000001.data
#cp ${datadir}/input_binaries_hires/pickup.0000103680.meta pickup.0000000001.meta
#cp ${datadir}/input_binaries_hires/pickup.0000103680.data pickup.0000000001.data
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

#### for HR adjoint
rm data.autodiff
cp ${rootdir}/input_adlo/data.autodiff .


#-- swap out data.ctrl and copy high-res adjustments
if [ ${iter} -lt 1 ]; then
  sed -i -e 's/'"doinitxx = .FALSE."'/'"doinitxx = .TRUE."'/g' data.ctrl
  sed -i -e 's/'"doInitXX = .FALSE."'/'"doInitXX = .TRUE."'/g' data.ctrl
  sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
  sed -i -e 's/'"doMainUnpack = .TRUE."'/'"doMainUnpack = .FALSE."'/g' data.ctrl
# starting from Argo assim
#  sed -i -e 's/'"doinitxx = .TRUE."'/'"doinitxx = .FALSE."'/g' data.ctrl
#  sed -i -e 's/'"doInitXX = .TRUE."'/'"doInitXX = .FALSE."'/g' data.ctrl
#  sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
#  sed -i -e 's/'"doMainUnpack = .TRUE."'/'"doMainUnpack = .FALSE."'/g' data.ctrl
#  cp ${datadir}/input_binaries_hires/xx_aqh.0000000006.data ./xx_aqh.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_aqh.effective.0000000006.data ./xx_aqh.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_atemp.0000000006.data ./xx_atemp.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_atemp.effective.0000000006.data ./xx_atemp.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_diffkr.0000000006.data ./xx_diffkr.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_diffkr.effective.0000000006.data ./xx_diffkr.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_lwdown.0000000006.data ./xx_lwdown.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_lwdown.effective.0000000006.data ./xx_lwdown.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_precip.0000000006.data ./xx_precip.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_precip.effective.0000000006.data ./xx_precip.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_salt.0000000006.data ./xx_salt.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_salt.effective.0000000006.data ./xx_salt.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_swdown.0000000006.data ./xx_swdown.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_swdown.effective.0000000006.data ./xx_swdown.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_theta.0000000006.data ./xx_theta.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_theta.effective.0000000006.data ./xx_theta.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_uwind.0000000006.data ./xx_uwind.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_uwind.effective.0000000006.data ./xx_uwind.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_vwind.0000000006.data ./xx_vwind.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_vwind.effective.0000000006.data ./xx_vwind.effective.0000000000.data 
#  cp ${datadir}/input_binaries_hires/xx_aqh.0000000006.meta ./xx_aqh.0000000000.meta
#  cp ${datadir}/input_binaries_hires/xx_aqh.effective.0000000006.meta ./xx_aqh.effective.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_atemp.0000000006.meta ./xx_atemp.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_atemp.effective.0000000006.meta ./xx_atemp.effective.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_diffkr.0000000006.meta ./xx_diffkr.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_diffkr.effective.0000000006.meta ./xx_diffkr.effective.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_lwdown.0000000006.meta ./xx_lwdown.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_lwdown.effective.0000000006.meta ./xx_lwdown.effective.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_precip.0000000006.meta ./xx_precip.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_precip.effective.0000000006.meta ./xx_precip.effective.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_salt.0000000006.meta ./xx_salt.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_salt.effective.0000000006.meta ./xx_salt.effective.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_swdown.0000000006.meta ./xx_swdown.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_swdown.effective.0000000006.meta ./xx_swdown.effective.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_theta.0000000006.meta ./xx_theta.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_theta.effective.0000000006.meta ./xx_theta.effective.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_uwind.0000000006.meta ./xx_uwind.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_uwind.effective.0000000006.meta ./xx_uwind.effective.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_vwind.0000000006.meta ./xx_vwind.0000000000.meta 
#  cp ${datadir}/input_binaries_hires/xx_vwind.effective.0000000006.meta ./xx_vwind.effective.0000000000.meta 
else
#  mkdir adxxfiles/
  cp ${dirrun_pup}/xx_hires/*000000${optimext}.data .
  cp -f ${datadir}/input_binaries_hires/xx_theta.0000000000.meta ./xx_theta.000000${optimext}.meta
  cp -f ${datadir}/input_binaries_hires/xx_atemp.0000000000.meta ./xx_atemp.000000${optimext}.meta
  cp $hr_dir_iter0/xx_theta.0000000000.meta ./xx_theta.000000${optimext}.meta
  cp $hr_dir_iter0/xx_salt.0000000000.meta ./xx_salt.000000${optimext}.meta
  cp $hr_dir_iter0/xx_atemp.0000000000.meta ./xx_atemp.000000${optimext}.meta
  cp $hr_dir_iter0/xx_precip.0000000000.meta ./xx_precip.000000${optimext}.meta
  cp $hr_dir_iter0/xx_swdown.0000000000.meta ./xx_swdown.000000${optimext}.meta
  cp $hr_dir_iter0/xx_lwdown.0000000000.meta ./xx_lwdown.000000${optimext}.meta
  cp $hr_dir_iter0/xx_uwind.0000000000.meta ./xx_uwind.000000${optimext}.meta
  cp $hr_dir_iter0/xx_vwind.0000000000.meta ./xx_vwind.000000${optimext}.meta
  cp $hr_dir_iter0/xx_aqh.0000000000.meta ./xx_aqh.000000${optimext}.meta
  cp $hr_dir_iter0/xx_diffkr.0000000000.meta ./xx_diffkr.000000${optimext}.meta
  sed -i -e 's/'"doinitxx = .TRUE."'/'"doinitxx = .FALSE."'/g' data.ctrl
  sed -i -e 's/'"doInitXX = .TRUE."'/'"doInitXX = .FALSE."'/g' data.ctrl
  sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
  sed -i -e 's/'"doMainUnpack = .FALSE."'/'"doMainUnpack = .FALSE."'/g' data.ctrl
  sed -i -e 's/'"doMainUnpack = .false."'/'"doMainUnpack = .false."'/g' data.ctrl

\rm data.optim
cat > data.optim <<EOF
 &OPTIM
 optimcycle=${iter},
 /
EOF

fi

#--- 7. executable --------
cp -p $builddir/mitgcmuv_ad .


set -x
date > run.MITGCM.timing
mpiexec -np ${nprocs_hr} ./mitgcmuv_ad > stdout
date >> run.MITGCM.timing


date > run.MITGCM.timing_adj
mpiexec -np ${nprocs_hr} ./mitgcmuv_ad > stdout_adj
date >> run.MITGCM.timing_adj


cd ..


