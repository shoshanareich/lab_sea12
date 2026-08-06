#!/bin/bash
#PBS -q long 
#PBS -l select=3:ncpus=40:model=sky_ele
#PBS -l walltime=48:00:00
#PBS -o argoMG.out
#PBS -e argoMG.err
#PBS -W umask=33
#PBS -m bea

#--- 0.load modules ------
ulimit -s unlimited
#limit stacksize unlimited
#module purge
#
#module load comp-intel/2018.3.222
#module load szip/2.1.1
#module load mpi-hpe/mpt
#module load hdf4/4.2.12
#module load hdf5/1.8.18_mpt
#module load netcdf/4.4.1.1_mpt
module use -a /swbuild/analytix/tools/modulefiles
module load miniconda3/v4

umask 022

ulimit -s hard
ulimit -u hard

#---- 1.set variables ------
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/lib
export FORT_BUFFERED=1
export MPI_IB_RAILS=2
export MPI_BUFS_PER_PROC=128
#with the current MPI, the following line doesn't matter.
#export MPI_BUFS_PER_HOST=512
export MPI_DISPLAY_SETTINGS


#---- set variables ------
# note: for nprocs, take ntiles - length(blanklist)
nprocs_hr=103
nprocs_lr=16

iter=0
itermax=10
costfactor=0.95

jobfile=run_argo_optimization_pfe.sh

#--- set dir ------------
rootdir=/nobackup/sreich/MITgcm_c69j/lab_sea12/
scratchdir=/nobackup/sreich/labsea_runs/assim_argo_MG/

builddir_hi=${rootdir}/build_adhi_2lev_seaice_update_mpi
builddir_lo=${rootdir}/build_adlo_2lev_seaice_update_mpi
optimdir=${scratchdir}/OPTIM

# --- optim ---

while [ ! ${iter} -gt $itermax ]; do

  ext2=$(printf "%04d" $iter)

# --- 1.high-res forward run ---  
  
  workdir_hi=${scratchdir}/run_adhi_it${ext2}

  if [ ! -d $workdir_hi ]; then
    mkdir -p $workdir_hi;
  fi
  
  cd $workdir_hi
  
  # cp binaries into workdir_hi
  # change data.ctrl
  # cp xx_[ctrl] adjustments if iter > 0
  ${rootdir}/link_hires_pfe.sh $iter $ext2 $scratchdir $builddir_hi 

  #---  run  --------
  \rm -f mitgcmuv*
  cp -f ${builddir_hi}/mitgcmuv_ad ./
  cp -f ${builddir_hi}/Makefile ./
  
  set -x
  date > run.MITGCM.timing
  mpiexec -np ${nprocs_hr} ./mitgcmuv_ad > stdout
  date >> run.MITGCM.timing
  cd ..

# --- 2.low-res adjoint run ---
  
  workdir_lo=${scratchdir}/run_adlo_it${ext2}

  if [ ! -d $workdir_lo ]; then
    mkdir -p $workdir_lo;
  fi

  source activate py38
  # create low-res xx_[ctrl]
  python3 ${rootdir}/mappings/make_cost_cp_v2.py "$ext2" "" "$scratchdir" "False"
  ## profiles retiling 
  python3 ${rootdir}/mappings/make_obsfit_lr_tiles.py "$ext2" "" "$scratchdir" "swot_obsfit_cycles_9thru11_labsea_L3v3" 
  python3 ${rootdir}/mappings/make_obsfit_lr_tiles.py "$ext2" "" "$scratchdir" "rads_20240101_20240308" 
#  python3 ${rootdir}/mappings/make_prof_lr_tiles.py "$ext2" "" "$scratchdir" "swot_obsfit_cycles_9thru11_labsea_L3v3_PROFILES"
  python3 ${rootdir}/mappings/make_prof_lr_tiles.py "$ext2" "" "$scratchdir" "ARGO_WO_2024_PFL_D_labsea_splitcost"

  cd $workdir_lo

  # cp binaries into workdir_lo
  # cp ONLINE low-res cost, misfit, barfiles, and xx_[ctrl]
  # create data.optim
  ${rootdir}/link_lores_pfe.sh $iter $workdir_hi $scratchdir $builddir_lo
  
  #---  run  --------
  \rm -f mitgcmuv*
  cp -f ${builddir_lo}/mitgcmuv_ad ./
  cp -f ${builddir_lo}/Makefile ./
  
  set -x
  date > run.MITGCM.timing
  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout
  date >> run.MITGCM.timing

#  sed -i 's/61/0/g' divided.ctrl
  sed -i 's/376/0/g' divided.ctrl
  date > run.MITGCM.timing
  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout
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

# --- 4. pack unpack -----------

  workdir_pup=${scratchdir}/run_adlo_packunpack

  if [ ! -d $workdir_pup ]; then
    mkdir -p $workdir_pup;
  fi

  cd $workdir_pup

  # cp ecco*ctrl*iter from optim 
  # change data.ctrl to unpack=TRUE
  # change data to run for 1 timestep
  # use lo build until figure out why we're not getting xx_[].effective
  # but it's only running the first chunk of the divided adjoint before automatically stopping 
  ${rootdir}/link_packunpack_pfe.sh $ext2 $workdir_pup $optimdir $builddir_lo $workdir_hi 

  #---  run  --------
  \rm -f mitgcmuv*
  cp -f ${builddir_lo}/mitgcmuv_ad ./
  cp -f ${builddir_lo}/Makefile ./
  
  set -x
  date > run.MITGCM.timing
  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout #2>&1 &

#  # Get the PID of the executable
#  EXEC_PID=$!
#  
#  # Monitor the log file in a loop using grep
#  while true; do
#      if grep -q "time_tsnumber" STDOUT.0000; then
#          echo "Pattern found, killing executable."
#          kill $EXEC_PID
#          break
#      fi
#      sleep 1  # Wait for 1 second before checking again
#  done
#  
#  # Optionally, wait for the executable to finish if not killed
#  wait $EXEC_PID

  date >> run.MITGCM.timing
  
  cd ..

# --- 5. interpolate adjustments to high-res -----------
python3 ${rootdir}/mappings/interp_xx_lores_to_hires_itX_v2.py $(printf "%04d" $((iter+1))) $workdir_pup  #"000$((iter+1))" $interp_type

#  let iter=6
  let iter=iter+1
done

echo "DONE"
