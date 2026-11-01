#!/bin/bash -x
#PBS -q long
#PBS -l select=3:ncpus=40:mpiprocs=40:model=sky_ele
#PBS -l walltime=24:00:00
#PBS -N swot10d_MGsLx1
#PBS -o /nobackup/sreich/MITgcm_c69j/lab_sea12/swot10d_MG_sL_x1_pfe.out
#PBS -e /nobackup/sreich/MITgcm_c69j/lab_sea12/swot10d_MG_sL_x1_pfe.err
#PBS -W umask=33
#PBS -m bea

# ============================================================================
# Pleiades repeat of submit_MG_10days_sigmaL.sh (sverdrup, x1 visc), to check
# that the two machines give the same optimization.  The science -- namelist
# seds, mapping arguments, optim settings, loop body -- is copied verbatim;
# only scheduler, paths, modules, conda and the MPI launcher differ.
#
# Submit directly:   qsub submit_MG_10days_sigmaL_pfe.sh
#
# Reference run on sverdrup (2026-10-09, 11 iterations, ~42 min each):
#   /scratch/shoshi/labsea_MG_12/assim_swot_MG_10days_sigmaL
#   costfinal it0000 1281566.04785003   it0010 146305.246026161
#
# Sizing: sky_ele is 40 cores/node, 103 HR ranks -> 3 nodes.  mpiprocs=40 is
# needed so MPT gets 120 rank slots, not 3.
#
# Inputs must be byte-identical to sverdrup's or the comparison means nothing;
# the md5 checks below enforce that (they are the sverdrup md5s).  Files to
# copy over first, if not already there:
#   input_binaries_hires/pickup{,_seaice}.0000103680.{data,meta}
#   input_binaries_hires/swot_obsfit_cycles_9thru11_labsea_L3v3_sigmaL.nc
#   input_binaries_lores/visc{d,z}_x1_bottom_topo.bin          (2026-10-08 17:52)
#   ${runsroot}/run_lo_spinup_viscx1_bottomtopo_fromrest/
#       pickup{,_seaice}.0000026298.{data,meta} and data       (the Oct 8 v2 rerun)
#   ${runsroot}/grid_{hires_cleanbathy,lores,lores_cleanbathy,lores_cleanbathy_v2}/
#   ${scratchdir}/OPTIM/optim.x
# Copy the LR spinup rather than rerunning it on pfe: a pfe spinup would differ
# at round-off and confound the cross-machine comparison.
# ============================================================================

# PBS runs this from a spool copy, so locate the checkout by absolute path.
export LABSEA_MACHINE=pfe
rootdir=/nobackup/sreich/MITgcm_c69j/lab_sea12
source ${rootdir}/machine_config.sh || exit 1
labsea_load_modules

iter=0
itermax=10
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

scratchdir=${runsroot}/assim_swot_MG_10days_sigmaL/
optimdir=${scratchdir}/OPTIM
workdir_pup=${scratchdir}/run_adlo_packunpack
builddir_hi=${rootdir}/build_adhi_2lev_seaice_update_mpi
builddir_lo=${rootdir}/build_adlo_2lev_seaice_update_mpi

spinup_lo=${runsroot}/run_lo_spinup_visc${visc}_bottomtopo_fromrest
pickup_lo=0000026298

for d in ${rootdir} ${scratchdir} ${optimdir} ${builddir_hi} ${builddir_lo} ${spinup_lo} \
         ${runsroot}/grid_hires_cleanbathy ${runsroot}/grid_lores \
         ${runsroot}/grid_lores_cleanbathy ${runsroot}/grid_lores_cleanbathy_v2; do
  if [ ! -d ${d} ]; then echo "FATAL: missing directory ${d}" >&2; exit 1; fi
done
for f in ${spinup_lo}/pickup.${pickup_lo}.data ${spinup_lo}/pickup_seaice.${pickup_lo}.data \
         ${rootdir}/input_binaries_lores/viscd_${visc}_bottom_topo.bin \
         ${rootdir}/input_binaries_lores/viscz_${visc}_bottom_topo.bin \
         ${rootdir}/input_binaries_hires/${obsfit_name}.nc \
         ${optimdir}/optim.x ${builddir_hi}/mitgcmuv_ad ${builddir_lo}/mitgcmuv_ad; do
  if [ ! -f ${f} ]; then echo "FATAL: missing ${f}" >&2; exit 1; fi
done

# sverdrup md5s of the run-specific inputs
md5_ok=true
while read sum f; do
  got=$(md5sum ${f} | cut -d' ' -f1)
  if [ "${got}" != "${sum}" ]; then
    echo "FATAL: ${f} differs from sverdrup (md5 ${got}, want ${sum})" >&2
    md5_ok=false
  fi
done <<EOF
c30a857aa0d69ad106ff18af4ef4abdb ${rootdir}/input_binaries_hires/pickup.0000103680.data
8da99d9c47fa13d60efe93e5356de545 ${rootdir}/input_binaries_hires/pickup_seaice.0000103680.data
78f36b92868bcb739e180e7788d7bf1b ${rootdir}/input_binaries_hires/${obsfit_name}.nc
361ac48a88a002a73f1495d258ea4eda ${rootdir}/input_binaries_lores/viscd_${visc}_bottom_topo.bin
03ff46579ef51dbd2907dab76141f602 ${rootdir}/input_binaries_lores/viscz_${visc}_bottom_topo.bin
dd56f375526b0ff4a369e1127243c2b6 ${spinup_lo}/pickup.${pickup_lo}.data
12df4482a37e22fcea3204f05e870b47 ${spinup_lo}/pickup_seaice.${pickup_lo}.data
EOF
${md5_ok} || exit 1

# record which code this ran, for the comparison
git -C ${rootdir} log -1 --format='%H %ci' > ${scratchdir}/git_version.txt
git -C ${rootdir} status --short >> ${scratchdir}/git_version.txt

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
  labsea_run ${nprocs_hr} ./mitgcmuv_ad > stdout_fwd
  date >> run.MITGCM.timing

  if [ ! -f ${workdir_hi}/costfinal ]; then
    echo "FATAL: iter ${iter} HR forward produced no costfinal" >&2; exit 1
  fi

# --- 2a. HR -> LR mappings ------------------------------------------------
  mkdir -p $workdir_lo/diags
  labsea_conda_activate py38
  cd ${rootdir}/mappings || exit 1
  python3 make_cost_cp_v2.py      "$ext2" "" "$scratchdir" "False" \
          $mp_factor $mp_nx $mp_ny $mp_maskdir $mp_restrict $mp_hrgrid
  python3 make_obsfit_lr_tiles.py "$ext2" "" "$scratchdir" "$obsfit_name"
  python3 make_obsfit_lr_tiles.py "$ext2" "" "$scratchdir" "rads_20240101_20240308"
  python3 make_prof_lr_tiles.py   "$ext2" "" "$scratchdir" "ARGO_WO_2024_PFL_D_labsea_splitcost"
  labsea_conda_deactivate

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
  labsea_run ${nprocs_lr} ./mitgcmuv_ad > stdout_adj1   # creates "377  376"
  sed -i 's/376/0/g' divided.ctrl                       # -> "377  0"
  labsea_run ${nprocs_lr} ./mitgcmuv_ad > stdout_adj2   # remaining levels
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
  labsea_run ${nprocs_lr} ./mitgcmuv_ad > stdout_pup

# --- 5. interpolate the LR adjustment up to HR ----------------------------
  labsea_conda_activate py38
  cd ${rootdir}/mappings || exit 1
  python3 interp_xx_lores_to_hires_itX_v2.py $(printf "%04d" $((iter+1))) $workdir_pup \
          $mp_factor $mp_nx $mp_ny $mp_grid $mp_method $mp_hrgrid
  labsea_conda_deactivate

  let iter=iter+1
done

echo "DONE"
