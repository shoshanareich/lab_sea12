#!/bin/bash -x
#SBATCH -J swot10d_MG_sL_x1_it11
#SBATCH -o swot10d_MG_sL_x1_it11.%j.out
#SBATCH -e swot10d_MG_sL_x1_it11.%j.err
#SBATCH -t 02:00:00
#SBATCH -N 6
#SBATCH -n 103
#SBATCH --mail-user=sreich@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# ============================================================================
# 10-day SWOT-only optimization.  HR forward (1/12) + 1/3 adjoint, multi-grid,
# injection prolongation + wet-area block-mean restriction (as in
# submit_MG_10days_newinterp.sh).
#
# Observations: swot_obsfit_cycles_9thru11_labsea_L3v3_sigmaL (representation
# error with L = dy ~ 37 km).  The name is passed to make_obsfit_lr_tiles.py,
# which globs the HR equi files on it -- the old name would silently re-tile
# nothing (or a stale set).
#
# LR adjoint viscosity: visc{d,z}_${visc}_bottom_topo.bin with the namelist of
# the matching spinup (viscAhGrid=0, viscA4GridMax=0.5, no viscAh), checked at
# run time against that spinup's data file.
#
# Pickups:
#   HR  run_hi_1yr_jra3q_B/pickup.0000103680                 360 d  (link_hires.sh)
#   LR  run_lo_spinup_visc${visc}_bottomtopo_fromrest/pickup.0000026298  365.25 d
# The 5.25-day HR/LR spinup-length mismatch is accepted deliberately.
#
# Sibling scripts differ only in visc= / scratchdir= / job name:
#   submit_MG_10days_sigmaL.sh          x1     (needs submit_LR_spinup_viscx1.sh first)
#   submit_MG_10days_sigmaL_viscx0.25.sh, _viscx0.1, _viscx0.01, _viscx0.001
# ============================================================================

module load intel impi netcdf netcdf-fortran phdf5
ulimit -s unlimited

iter=11              # one more x1 iteration, in place (warm start from OPWARM.opt0011)
itermax=11
nprocs_hr=103
nprocs_lr=16

visc=x1
nts_hr=2880          # 10 days at deltaT=300
nts_lr=720           # 10 days at deltaT=1200
msfac=4              # 1/12 -> 1/3
obsfit_name=swot_obsfit_cycles_9thru11_labsea_L3v3_sigmaL

# resolution + operator arguments handed to mappings/*.py
mp_factor=4
mp_nx=120
mp_ny=96
mp_grid=grid_lores_cleanbathy_v2
mp_hrgrid=grid_hires_cleanbathy
mp_method=injection   # piecewise-constant prolongation (LR -> HR)
mp_maskdir=grid_lores_cleanbathy_v2
mp_restrict=area      # wet-area-weighted block mean (HR -> LR); exact adjoint of injection

rootdir=/home/shoshi/MITgcm_c69j/lab_sea12/
scratchdir=/scratch/shoshi/labsea_MG_12/assim_swot_MG_10days_sigmaL/
optimdir=${scratchdir}/OPTIM
workdir_pup=${scratchdir}/run_adlo_packunpack
builddir_hi=${rootdir}/build_adhi_2lev_seaice_update_mpi
builddir_lo=${rootdir}/build_adlo_2lev_seaice_update_mpi

spinup_lo=/scratch/shoshi/labsea_MG_12/run_lo_spinup_visc${visc}_bottomtopo_fromrest
pickup_lo=0000026298

for d in ${rootdir} ${scratchdir} ${optimdir} ${builddir_hi} ${builddir_lo} ${spinup_lo}; do
  if [ ! -d ${d} ]; then echo "FATAL: missing directory ${d}" >&2; exit 1; fi
done
for f in ${spinup_lo}/pickup.${pickup_lo}.data ${spinup_lo}/pickup_seaice.${pickup_lo}.data \
         ${rootdir}/input_binaries_lores/viscd_${visc}_bottom_topo.bin \
         ${rootdir}/input_binaries_lores/viscz_${visc}_bottom_topo.bin \
         ${rootdir}/input_binaries_hires/${obsfit_name}.nc; do
  if [ ! -f ${f} ]; then echo "FATAL: missing ${f}" >&2; exit 1; fi
done

# SWOT-only cost + 10-day window + tape coarsening, applied in $PWD.
#   $1 = nTimeSteps for this resolution   $2 = MSfacX/MSfacY
set_namelists() {
  local nts=$1 msf=$2
  sed -i -e "s/^ *obsfitFiles(1) *=.*/ obsfitFiles(1) = '${obsfit_name}',/g" data.obsfit
  sed -i -e 's/^ *mult_obsfit(1) *=.*/ mult_obsfit(1) = 1.0,/g' data.obsfit
  sed -i -e 's/^ *mult_obsfit(2) *=.*/ mult_obsfit(2) = 0.0,/g' data.obsfit
  sed -i -e 's/^ *mult_profiles(1,1) *=.*/ mult_profiles(1,1)   = 0.,/g' data.profiles
  sed -i -e 's/^ *mult_profiles(1,2) *=.*/ mult_profiles(1,2)   = 0.,/g' data.profiles
  sed -i -e 's/^ *mult_profiles_mean(1) *=.*/ mult_profiles_mean(1)  = 0.,/g' data.profiles
  sed -i -e "s/^ *nTimeSteps *=.*/ nTimeSteps=${nts},/g" data
  sed -i -e "s/^ *nIter0 *=.*/ nIter0=1,/g" data   # pickup*.0000000001
  sed -i -e "s/^ *MSfacX *=.*/ MSfacX = ${msf},/g" data.autodiff
  sed -i -e "s/^ *MSfacY *=.*/ MSfacY = ${msf},/g" data.autodiff
  if ! grep -q "^ *obsfitFiles(1) *= *'${obsfit_name}'" data.obsfit; then
    echo "FATAL: data.obsfit in $PWD does not point at ${obsfit_name}" >&2; exit 1
  fi
}

# LR only: viscosity namelist + spinup pickup, applied in $PWD after link_*.sh.
# The physics lines must match the spinup the pickup came from.
phys_lines() { grep -E '^ *(visc|diffK|tempAdv|saltAdv|no_slip|implicitVisc)' $1; }
set_lr_visc_pickup() {
  sed -i -e 's/^ *viscAh *=.*/# viscAh=1.E0,/' \
         -e 's/^ *viscAhGrid *=.*/ viscAhGrid = 0.E-4,/' \
         -e "s/^ *viscAHDfile *=.*/ viscAHDfile    ='viscd_${visc}_bottom_topo.bin',/" \
         -e "s/^ *viscAHZfile *=.*/ viscAHZfile    ='viscz_${visc}_bottom_topo.bin',/" data
  grep -q '^ *viscA4GridMax *=' data || sed -i '/^ *viscC4Leithd *=/a\ viscA4GridMax=0.5,' data
  if ! diff <(phys_lines data) <(phys_lines ${spinup_lo}/data) >&2; then
    echo "FATAL: LR physics namelist in $PWD differs from ${spinup_lo}/data" >&2; exit 1
  fi
  # link_*.sh only link visc[dz]_bottom_topo*.bin, which misses the x-series names
  for c in d z; do
    ln -sf ${rootdir}/input_binaries_lores/visc${c}_${visc}_bottom_topo.bin .
    if [ ! -e visc${c}_${visc}_bottom_topo.bin ]; then
      echo "FATAL: visc${c}_${visc}_bottom_topo.bin not staged in $PWD" >&2; exit 1
    fi
  done
  # link_lores.sh / link_packunpack.sh stage run_lo_1yr_jra3q_B's 25920; replace
  for p in pickup pickup_seaice; do
    for s in meta data; do
      cp -f ${spinup_lo}/${p}.${pickup_lo}.${s} ./${p}.0000000001.${s}
    done
  done
}

while [ ! ${iter} -gt $itermax ]; do

  ext2=$(printf "%04d" $iter)
  workdir_hi=${scratchdir}/run_adhi_it${ext2}
  workdir_lo=${scratchdir}/run_adlo_it${ext2}

# --- 1. HR forward (writes coarsened tapes, cost, misfits) -----------------
  mkdir -p $workdir_hi
  cd $workdir_hi || exit 1

  # real iter here: link_hires.sh's iter>0 branch pulls xx_hires/ from the
  # pack/unpack dir, which is exactly the multi-grid path we want.
  ${rootdir}/link_hires.sh $iter $ext2 $scratchdir $builddir_hi
  set_namelists $nts_hr $msfac

  \rm -f costfinal divided.ctrl
  \rm -f mitgcmuv*
  cp -f ${builddir_hi}/mitgcmuv_ad ./
  date > run.MITGCM.timing
  mpiexec -np ${nprocs_hr} ./mitgcmuv_ad > stdout_fwd
  date >> run.MITGCM.timing

  if [ ! -f ${workdir_hi}/costfinal ]; then
    echo "FATAL: iter ${iter} HR forward produced no costfinal" >&2; exit 1
  fi

# --- 2a. HR -> LR mappings ------------------------------------------------
  mkdir -p $workdir_lo/diags
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate py38
  cd ${rootdir}/mappings || exit 1
  python3 make_cost_cp_v2.py      "$ext2" "" "$scratchdir" "False" \
          $mp_factor $mp_nx $mp_ny $mp_maskdir $mp_restrict $mp_hrgrid
  python3 make_obsfit_lr_tiles.py "$ext2" "" "$scratchdir" "$obsfit_name"
  python3 make_obsfit_lr_tiles.py "$ext2" "" "$scratchdir" "rads_20240101_20240308"
  python3 make_prof_lr_tiles.py   "$ext2" "" "$scratchdir" "ARGO_WO_2024_PFL_D_labsea_splitcost"
  conda deactivate

  if ! ls ${workdir_hi}/obs_LR_equi/${obsfit_name}.*.equi.nc > /dev/null 2>&1; then
    echo "FATAL: iter ${iter}: no LR-tiled ${obsfit_name} equi files" >&2; exit 1
  fi

# --- 2b. LR adjoint -------------------------------------------------------
  cd $workdir_lo || exit 1
  ${rootdir}/link_lores.sh $iter $workdir_hi $scratchdir $builddir_lo
  set_namelists $nts_lr 1
  set_lr_visc_pickup

  \rm -f divided.ctrl
  \rm -f mitgcmuv*
  cp -f ${builddir_lo}/mitgcmuv_ad ./
  date > run.MITGCM.timing
  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout_adj1   # creates "377  376"
  sed -i 's/376/0/g' divided.ctrl                        # -> "377  0"
  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout_adj2   # remaining levels
  date >> run.MITGCM.timing

  if [ ! -f ${workdir_lo}/ecco_cost_MIT_CE_000.opt${ext2} ]; then
    echo "FATAL: iter ${iter} produced no gradient file; stopping" >&2; exit 1
  fi

# --- 3. OPTIM (at LR) -----------------------------------------------------
  cd ${optimdir} || exit 1
  cp ${workdir_lo}/ecco_cost_MIT_CE_000.opt${ext2} ${optimdir}
  cp ${workdir_lo}/ecco_ctrl_MIT_CE_000.opt${ext2} ${optimdir}
  cp ${workdir_lo}/data.ctrl ${optimdir}
  cp -f ${workdir_lo}/costfinal ${optimdir}
  echo "iter = $iter"

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
    iprint=10,
    nupdate=4,
    /
    &M1QN3
    coldstart = ${coldstart},
    /
EOF
  ./optim.x > opt_it${ext2}.txt
  cp data.optim data.optim_${ext2}

# --- 4. pack/unpack the new LR control ------------------------------------
  mkdir -p $workdir_pup
  cd $workdir_pup || exit 1
  ${rootdir}/link_packunpack.sh $ext2 $workdir_pup $optimdir $builddir_lo $workdir_hi
  set_namelists $nts_lr 1
  set_lr_visc_pickup

  \rm -f mitgcmuv*
  cp -f ${builddir_lo}/mitgcmuv_ad ./
  mpiexec -np ${nprocs_lr} ./mitgcmuv_ad > stdout_pup

# --- 5. interpolate the LR adjustment up to HR ----------------------------
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate py38
  cd ${rootdir}/mappings || exit 1
  python3 interp_xx_lores_to_hires_itX_v2.py $(printf "%04d" $((iter+1))) $workdir_pup \
          $mp_factor $mp_nx $mp_ny $mp_grid $mp_method $mp_hrgrid
  conda deactivate

  let iter=iter+1
done

echo "DONE"
