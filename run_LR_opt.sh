#!/bin/bash -x
# ---------------------------------------------------------------------------
# run_LR_opt.sh
#
# Portable across sverdrup / stampede3 / pfe: every machine-specific path,
# module load, MPI launcher and conda activation comes from machine_config.sh.
# Override the detected machine with:
#     export LABSEA_MACHINE=sverdrup|stampede3|pfe
#
# Submit with the wrapper for your machine:
#     ./submit.sh run_LR_opt.sh          # detects machine, picks sbatch or qsub
# or run directly inside an existing allocation:
#     ./run_LR_opt.sh
# ---------------------------------------------------------------------------

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/machine_config.sh"
labsea_load_modules

#---- set variables ------
# nprocs_hr / nprocs_lr default to 103 / 16 in machine_config.sh.
# They follow SIZE.h + data.exch2 of the build, not the machine:
#   nprocs = (number of tiles) - (length of blankList)
# Override below only if you point builddir_* at a different tiling.

iter=6
itermax=20
costfactor=0.95

jobfile=run_optimization.bash

#--- set dir ------------
scratchdir=${runsroot}/assim_swot_LR

builddir_lo=${rootdir}/build_adlo_2lev_seaice_update_mpi
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
  cp ${rootdir}/input_adlo/* .
  ln -s ${rootdir}/input_binaries_lores/bathy_cleaned_96x120.bin .
  ln -s ${rootdir}/input_binaries_lores/LevCli_temp_120x96_linearv3_smooth.labsea1979 .
  ln -s ${rootdir}/input_binaries_lores/LevCli_salt_120x96_linearv3_smooth.labsea1979 .
  ln -s ${rootdir}/input_binaries_lores/viscd_bottom_topo*.bin .
  ln -s ${rootdir}/input_binaries_lores/viscz_bottom_topo*.bin .
  ln -s ${rootdir}/input_binaries_lores/diffkr_r4.bin .
  ln -s ${rootdir}/input_binaries_exf/* .
  ln -s ${rootdir}/input_binaries_hires/ones_64b.bin .
  ln -s ${rootdir}/input_binaries_hires/ARGO_WO_2024_PFL_D_labsea_splitcost.nc .
  ln -s ${rootdir}/input_binaries_hires/swot*.nc .
  ln -s ${rootdir}/input_binaries_hires/rads_20240101_20240308.nc .
  ln -s ${rootdir}/input_binaries_hires/rads*.nc .
  ln -s ${rootdir}/input_binaries_lores/rads_*_v2_2024* .
  ln -s ${rootdir}/input_weights_lores/*_jra3q_weights_Jan2024_64b_SMOOTHED_removeboundary.bin .
  ln -s ${rootdir}/input_binaries_lores/rads_j3_labsea_96x120_v2_2024 .
  ln -s ${rootdir}/input_binaries_lores/slaerr_03m.bin .
  ln -s ${rootdir}/input_weights_lores/*fromASTE_*.bin .
  cp ${rootdir}/input_binaries_lores/pickup_seaice.0000025920.meta pickup_seaice.0000000001.meta
  cp ${rootdir}/input_binaries_lores/pickup_seaice.0000025920.data pickup_seaice.0000000001.data
  cp ${rootdir}/input_binaries_lores/pickup.0000025920.meta pickup.0000000001.meta
  cp ${rootdir}/input_binaries_lores/pickup.0000025920.data pickup.0000000001.data
  
  mkdir jra3q
  #labsea_link_jra3q ./jra3q
  rm data.diagnostics
  cp ${rootdir}/input_adhi/data.diagnostics .
  
  ##--- 5. linking xx_ fields ------
  if [ ${iter} -lt 1 ]; then
    sed -i -e 's/'"doinitxx = .FALSE."'/'"doinitxx = .FALSE."'/g' data.ctrl
    sed -i -e 's/'"doInitXX = .FALSE."'/'"doInitXX = .FALSE."'/g' data.ctrl
    sed -i -e 's/'"doMainPack = .FALSE."'/'"doMainPack = .TRUE."'/g' data.ctrl
    sed -i -e 's/'"doMainUnpack = .TRUE."'/'"doMainUnpack = .FALSE."'/g' data.ctrl
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_aqh.0000000009.data ./xx_aqh.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_aqh.effective.0000000009.data ./xx_aqh.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_atemp.0000000009.data ./xx_atemp.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_atemp.effective.0000000009.data ./xx_atemp.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_diffkr.0000000009.data ./xx_diffkr.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_diffkr.effective.0000000009.data ./xx_diffkr.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_lwdown.0000000009.data ./xx_lwdown.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_lwdown.effective.0000000009.data ./xx_lwdown.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_precip.0000000009.data ./xx_precip.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_precip.effective.0000000009.data ./xx_precip.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_salt.0000000009.data ./xx_salt.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_salt.effective.0000000009.data ./xx_salt.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_swdown.0000000009.data ./xx_swdown.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_swdown.effective.0000000009.data ./xx_swdown.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_theta.0000000009.data ./xx_theta.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_theta.effective.0000000009.data ./xx_theta.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_uwind.0000000009.data ./xx_uwind.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_uwind.effective.0000000009.data ./xx_uwind.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_vwind.0000000009.data ./xx_vwind.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_vwind.effective.0000000009.data ./xx_vwind.effective.0000000000.data
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_aqh.0000000009.meta ./xx_aqh.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_aqh.effective.0000000009.meta ./xx_aqh.effective.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_atemp.0000000009.meta ./xx_atemp.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_atemp.effective.0000000009.meta ./xx_atemp.effective.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_diffkr.0000000009.meta ./xx_diffkr.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_diffkr.effective.0000000009.meta ./xx_diffkr.effective.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_lwdown.0000000009.meta ./xx_lwdown.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_lwdown.effective.0000000009.meta ./xx_lwdown.effective.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_precip.0000000009.meta ./xx_precip.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_precip.effective.0000000009.meta ./xx_precip.effective.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_salt.0000000009.meta ./xx_salt.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_salt.effective.0000000009.meta ./xx_salt.effective.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_swdown.0000000009.meta ./xx_swdown.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_swdown.effective.0000000009.meta ./xx_swdown.effective.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_theta.0000000009.meta ./xx_theta.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_theta.effective.0000000009.meta ./xx_theta.effective.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_uwind.0000000009.meta ./xx_uwind.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_uwind.effective.0000000009.meta ./xx_uwind.effective.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_vwind.0000000009.meta ./xx_vwind.0000000000.meta
    cp ${runsroot}/assim_argo_LR/run_lo_it0009/xx_vwind.effective.0000000009.meta ./xx_vwind.effective.0000000000.meta
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
  labsea_run ${nprocs_lr} ./mitgcmuv_ad > stdout # run forward
  labsea_run ${nprocs_lr} ./mitgcmuv_ad > stdout # run first chunk of DivA

  sed -i 's/376/0/g' divided.ctrl 
  labsea_run ${nprocs_lr} ./mitgcmuv_ad > stdout # run rest of DivA
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
