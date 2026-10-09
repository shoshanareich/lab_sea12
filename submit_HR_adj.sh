#!/bin/bash -x
#SBATCH -J labseaHRswot_adj
#SBATCH -o labseaHRswot_adj.%j.out
#SBATCH -e labseaHRswot_adj.%j.err
#SBATCH -t 48:00:00
#SBATCH -N 6
#SBATCH -n 103
#SBATCH --mail-user=sreich@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# ============================================================================
# Single-resolution 1/12 deg (HR) adjoint optimization -- SWOT ONLY.
#

#--- 0. modules / environment (sverdrup) ---
module load intel impi netcdf netcdf-fortran phdf5
ulimit -s unlimited

iter=0
itermax=3
nprocs_hr=103

jobfile=submit_HR_adj.sh

#--- 2.set dir ------------
rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
scratchdir=/scratch/shoshi/labsea_MG_12/assim_swot_HR/
optimdir=${scratchdir}/OPTIM

builddir=${rootdir}/build_adhi_2lev_seaice_update_mpi

# Fail loudly rather than looping through with empty paths.
for d in ${rootdir} ${scratchdir} ${builddir}; do
  if [ ! -d ${d} ]; then echo "FATAL: missing directory ${d}" >&2; exit 1; fi
done
if [ ! -x ${builddir}/mitgcmuv_ad ]; then
  echo "FATAL: no mitgcmuv_ad in ${builddir}" >&2; exit 1
fi

mkdir -p ${optimdir}

while [ ! ${iter} -gt $itermax ]; do

  ext2=$(printf "%04d" $iter)

# --- 1.high-res forward run ---

  workdir_hi=${scratchdir}/run_adhi_it${ext2}

  if [ ! -d $workdir_hi ]; then
    mkdir -p $workdir_hi;
  fi

  cd $workdir_hi || exit 1

  # link_hires.sh is called with iter=0 on purpose: its iter>0 branch is the
  # multi-grid path (xx_hires/ out of run_adlo_packunpack), which does not
  # exist here.  Its data.ctrl flags are overridden just below.
  ${rootdir}/link_hires.sh 0 $ext2 $scratchdir $builddir

  #--- namelist overrides ---------------------------------------------------
  if [ ${iter} -eq 0 ]; then
    sed -i -e 's/^ *doInitXX *=.*/ doInitXX = .TRUE.,/g' data.ctrl
  else
    sed -i -e 's/^ *doInitXX *=.*/ doInitXX = .FALSE.,/g' data.ctrl
  fi
  sed -i -e 's/^ *doMainPack *=.*/ doMainPack = .TRUE.,/g' data.ctrl
  sed -i -e 's/^ *doMainUnpack *=.*/ doMainUnpack = .TRUE.,/g' data.ctrl

  # data.obsfit: SWOT on.  input_adhi/data.obsfit ships with mult_obsfit(1)=0.0,
  # so this has to be set explicitly or the cost has no observations at all.
  sed -i -e 's/^ *mult_obsfit(1) *=.*/ mult_obsfit(1) = 1.0,/g' data.obsfit

  # data.profiles: Argo OFF for the SWOT-only experiment.
  sed -i -e 's/^ *mult_profiles(1,1) *=.*/ mult_profiles(1,1)   = 0.,/g' data.profiles
  sed -i -e 's/^ *mult_profiles(1,2) *=.*/ mult_profiles(1,2)   = 0.,/g' data.profiles

  # data.autodiff: MSfacX=MSfacY=1.  input_adhi ships 4,4 for the multi-grid
  # case, where the HR forward writes tapes coarsened x4 for the 1/3 deg
  # adjoint.  A single-resolution HR adjoint must read full-size tapes.
  sed -i -e 's/^ *MSfacX *=.*/ MSfacX = 1,/g' data.autodiff
  sed -i -e 's/^ *MSfacY *=.*/ MSfacY = 1,/g' data.autodiff

  #--- (re)set optimcycle ---------------------------------------------------
  \rm -f data.optim
  cat > data.optim <<EOF
 &OPTIM
 optimcycle=${iter},
 /
EOF

  # new control vector from the previous OPTIM (unpacked by the model)
  if [ ${iter} -gt 0 ]; then
    cp ${optimdir}/ecco_ctrl_MIT_CE_000.opt${ext2} .
  fi

  # Clear leftover DivA state.  A stale costfinal suppresses CTRL_UNPACK, and a
  # stale divided.ctrl makes the three-invocation sequence below do the wrong
  # number of levels.
  \rm -f costfinal divided.ctrl

  #---  run  --------
  \rm -f mitgcmuv*
  cp -f ${builddir}/mitgcmuv_ad ./
  cp -f ${builddir}/Makefile ./

  set -x
  date > run.MITGCM.timing
  mpiexec -np ${nprocs_hr} ./mitgcmuv_ad > stdout_fwd   # forward: cost + tapes
  date >> run.MITGCM.timing

# --- 2.high-res adjoint run ---
  date >> run.MITGCM.timing
  mpiexec -np ${nprocs_hr} ./mitgcmuv_ad > stdout_adj1  # creates "377  376"

#  sed -i 's/61/0/g' divided.ctrl
  sed -i 's/376/0/g' divided.ctrl                       # -> "377  0"
  mpiexec -np ${nprocs_hr} ./mitgcmuv_ad > stdout_adj2  # remaining DivA levels
  date >> run.MITGCM.timing

  if [ ! -f ${workdir_hi}/ecco_cost_MIT_CE_000.opt${ext2} ]; then
    echo "FATAL: iter ${iter} produced no gradient file; stopping" >&2; exit 1
  fi


# --- 3. OPTIM ----------------
  cd ${optimdir} || exit 1
  #bash reset.bash
  cp ${workdir_hi}/ecco_cost_MIT_CE_000.opt${ext2} ${optimdir}
  cp ${workdir_hi}/ecco_ctrl_MIT_CE_000.opt${ext2} ${optimdir}
  cp ${workdir_hi}/data.ctrl ${optimdir}
  cp -f ${workdir_hi}/costfinal ${optimdir}
  echo "iter = $iter"

  # Warm-start m1qn3 whenever the previous iteration left an OPWARM file, so
  # the quasi-Newton memory carries across outer iterations.  OPWARM.opt<iter>
  # is NOT deleted here for the same reason.
  if [ -f ${optimdir}/OPWARM.opt${ext2} ]; then
    coldstart=".FALSE."
  else
    coldstart=".TRUE."
  fi

  mv -f data.optim data.optim_bk
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
    coldstart = ${coldstart},
    /
EOF

  ./optim.x > opt_it${ext2}.txt
  cp data.optim data.optim_${ext2}


  let iter=iter+1
done

echo "DONE"
