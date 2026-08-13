#!/bin/bash -x
#SBATCH -J swot_LR6
#SBATCH -o swot_LR6.%j.out
#SBATCH -e swot_LR6.%j.err
#SBATCH -t 48:00:00
#SBATCH -N 2 
#SBATCH -n 48


#--- 0.load modules ------

#---- set variables ------
# note: for nprocs, take ntiles - length(blanklist)
nprocs_lr=48

iter=0
itermax=1
costfactor=0.95

jobfile=run_LR6_opt.sh

#--- set dir ------------
rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
scratchdir=/scratch/shoshi/labsea_MG_12/assim_swot_LR6/

builddir_lo=${rootdir}/build_adlo6_2lev_seaice_update_mpi
optimdir=${scratchdir}/OPTIM

# --- optim ---

while [ ! ${iter} -gt $itermax ]; do

  ext2=$(printf "%04d" $iter)

# --- low-res forward and adjoint run ---
  
  workdir_lo=${scratchdir}/run_lo_it${ext2}

  if [ ! -d $workdir_lo ]; then
    mkdir -p $workdir_lo;
  fi

  cd $workdir_lo

  mkdir -p ./diags/     

  #--- 6. NAMELISTS ---------
  #ln -s ${rootdir}/input_cal/* .
  cp ${rootdir}/input_adlo6/* .
  ln -s ${rootdir}/input_binaries_lores6/bathy_cleaned_192x240.bin .
  ln -s ${rootdir}/input_binaries_lores6/LevCli_temp_240x192_v1_binavg.labsea1979 .
  ln -s ${rootdir}/input_binaries_lores6/LevCli_salt_240x192_v1_binavg.labsea1979 .
  ln -s ${rootdir}/input_binaries_lores6/visc*.bin .
  ln -s ${rootdir}/input_binaries_lores6/visc*.bin .
  ln -s ${rootdir}/input_binaries_lores6/diffkr_r4_240x192.bin .
  ln -s ${rootdir}/input_binaries_lores6/KapRedi_from_aste540* .
  ln -s ${rootdir}/input_binaries_lores6/KapGM_from_aste540* .
  ln -s ${rootdir}/input_binaries_exf/* .
  ln -s ${rootdir}/input_binaries_hires/ones_64b.bin .
  ln -s ${rootdir}/input_binaries_hires/ARGO_WO_2024_PFL_D_labsea_splitcost.nc .
  ln -s ${rootdir}/input_binaries_hires/swot*.nc .
  ln -s ${rootdir}/input_binaries_hires/rads_20240101_20240308.nc .
  ln -s ${rootdir}/input_binaries_hires/rads*.nc .
  ln -s ${rootdir}/input_binaries_lores6/rads_*_2024* .
  ln -s ${rootdir}/input_binaries_lores6/slaerr_3cm_192x240.bin .
  ln -s ${rootdir}/input_weights_lores6/* .
  cp ${rootdir}/input_binaries_lores6/pickup_seaice.0000051840.meta pickup_seaice.0000000001.meta
  cp ${rootdir}/input_binaries_lores6/pickup_seaice.0000051840.data pickup_seaice.0000000001.data
  cp ${rootdir}/input_binaries_lores6/pickup.0000051840.meta pickup.0000000001.meta
  cp ${rootdir}/input_binaries_lores6/pickup.0000051840.data pickup.0000000001.data
  
  mkdir jra3q
  ln -s /scratch/shared/jra3q/jra3q_*_2024 ./jra3q/
  #ln -s /scratch/08382/shoshi/jra3q/jra3q_*_2024 ./jra3q/
  rm data.diagnostics
  cp ${rootdir}/input_adhi/data.diagnostics .
#  cp /scratch/shoshi/labsea_MG_12/test_lores/run_spinup_1yr_C/data.diagnostics .

  ##--- 5. linking xx_ fields ------
  if [ ${iter} -lt 1 ]; then
    sed -i -e 's/'"doinitxx = .FALSE."'/'"doinitxx = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doInitXX = .FALSE."'/'"doInitXX = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doMainUnpack = .TRUE."'/'"doMainUnpack = .FALSE."'/g' data.ctrl
  else
    sed -i -e 's/'"doinitxx = .TRUE."'/'"doinitxx = .FALSE."'/g' data.ctrl
    sed -i -e 's/'"doInitXX = .TRUE."'/'"doInitXX = .FALSE."'/g' data.ctrl
    sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doMainUnpack = .FALSE."'/'"doMainUnpack = .TRUE."'/g' data.ctrl
    cp ${optimdir}/ecco_ctrl_MIT_CE_000.opt${ext2} .
  fi
  #--- 10. (re)set optimcycle --------------------
  
  \rm data.optim
  cat > data.optim <<EOF
   &OPTIM
   optimcycle=${iter},
   /
EOF

  #---  run forward  --------
  cp -f ${builddir_lo}/mitgcmuv_ad ./
  cp -f ${builddir_lo}/Makefile ./
  
  set -x
  date > run.MITGCM.timing
  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout # run forward
  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout # run first chunk of DivA

  sed -i 's/376/0/g' divided.ctrl 
  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout # run rest of DivA
  date >> run.MITGCM.timing
  cd ..


# --- 3. OPTIM ----------------
  cd ${optimdir}
  #bash reset.bash
  cp ${workdir_lo}/ecco_cost_MIT_CE_000.opt${ext2} ${optimdir} 
  cp ${workdir_lo}/ecco_ctrl_MIT_CE_000.opt${ext2} ${optimdir} 
  cp ${workdir_lo}/data.ctrl ${optimdir} 
  cp -f ${workdir_lo}/costfinal ${optimdir}
#  cost=$(grep fc costfunction${ext2}  | sed 's/D/E/g' | awk '{printf "%14.12e", $3}')
#  costf=$(grep fc costfunction${ext2} | sed 's/D/E/g' | awk '{printf "%0.14f", $3}')
  echo "iter = $iter"
#  echo "cost = $cost"
  
#  costupdate=$(echo $costf*$costfactor | bc)
#  costnew=`echo $costupdate|awk '{printf "%14.12e\n", $costupdate}'`
#  echo "costnew = $costnew"

  mv data.optim data.optim_bk
  cat > data.optim <<EOF
    &OPTIM
    optimcycle=${iter},
    numiter=100,
    nfunc=9,
    dfminFrac=0.05,
#    fmin=${costnew},
    iprint=10,
    nupdate=4,
    /
    &M1QN3
    coldstart = .TRUE.,
    /
EOF

  \rm OP*
  ./optim.x > opt_it${ext2}.txt

  dir_iter=${workdir_lo}/optim/+it${ext2}
  mkdir -p $dir_iter
  cp data.optim $dir_iter
  cp data.optim data.optim_${ext2}
  cp ecco_ctrl_MIT_CE_000.opt${ext2} $dir_iter

  cd ..

  let iter=iter+1
done

echo "DONE"
