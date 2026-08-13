#!/bin/bash
# ============================================================================
# machine_config.sh -- shared per-machine configuration for lab_sea12
#
# Source this near the top of any run_*.sh or link_*.sh:
#
#     source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/machine_config.sh"
#
# Auto-detects sverdrup / stampede3 / pfe.  To override (e.g. on an odd
# compute node, or to test another machine's paths):
#
#     export LABSEA_MACHINE=sverdrup|stampede3|pfe
#
# Exports:
#   LABSEA_MACHINE  resolved machine name
#   rootdir         lab_sea12 repo root -- code_*, input_ad*, build_*, mappings/,
#                   and input_binaries_* / input_weights_* (symlinked on sverdrup)
#   runsroot        scratch root holding run_* dirs and grid_hires/grid_lores
#   mitgcmutils     path to MITgcm/utils/python/MITgcmutils
#   optfile         genmake2 -optfile argument
#   scheduler       slurm | pbs
#   mpi_launch      "mpiexec -np" or "ibrun -n"
#   nprocs_hr       HR rank count (ntiles - blanks; 32x32 tiling => 103)
#   nprocs_lr       LR rank count (5x4 tiles - 4 blanks => 16)
#
# Functions:
#   labsea_run <nprocs> <exe> ...   launch under the right MPI launcher
#   labsea_load_modules             module purge/load for this machine
#   labsea_conda_activate [env]     default env py38
#   labsea_conda_deactivate
#   labsea_link_jra3q [destdir]     default ./jra3q
# ============================================================================

# ---------------------------------------------------------------- 1. detect --
if [ -n "${LABSEA_MACHINE:-}" ]; then
    _labsea_how="LABSEA_MACHINE override"   # explicit override wins
elif [ -d /nobackup/sreich ]; then
    LABSEA_MACHINE=pfe;       _labsea_how="found /nobackup/sreich"
elif [ -n "${TACC_SYSTEM:-}" ] || [ -d /work2/08382/shoshi ]; then
    LABSEA_MACHINE=stampede3; _labsea_how="TACC_SYSTEM=${TACC_SYSTEM:-unset} / found /work2/08382/shoshi"
elif [ -d /scratch/shoshi ] || [[ "$(hostname -f)" == sverdrup* ]]; then
    LABSEA_MACHINE=sverdrup;  _labsea_how="found /scratch/shoshi or hostname $(hostname -s)"
else
    echo "machine_config.sh: ERROR -- cannot identify this machine." >&2
    echo "  Set LABSEA_MACHINE=sverdrup|stampede3|pfe and re-run." >&2
    return 1 2>/dev/null || exit 1
fi
export LABSEA_MACHINE

# Announce which machine we resolved to, and why. This is the single most
# useful line when a run stages the wrong files or lands in the wrong scratch.
# Silence it with: export LABSEA_QUIET=1
if [ -z "${LABSEA_QUIET:-}" ]; then
    echo "[machine_config] LABSEA_MACHINE=${LABSEA_MACHINE}  (${_labsea_how})"
fi
unset _labsea_how

# --------------------------------------------------- 2. per-machine settings --
case "${LABSEA_MACHINE}" in

  sverdrup)
    rootdir=/home/shoshi/MITgcm_c69j/lab_sea12
    runsroot=/scratch/shoshi/labsea_MG_12
    mitgcmutils=/home/shoshi/MITgcm_c69j/MITgcm/utils/python/MITgcmutils
    optfile=/home/shoshi/computing/optfiles/linux_amd64_ifort+mpi_sverdrup
    scheduler=slurm
    mpi_launch="mpiexec -np"
    jra3q_src=( '/scratch/shared/jra3q/jra3q_*_2024' )
    ;;

  stampede3)
    rootdir=/work2/08382/shoshi/stampede3/MITgcm_c69j/lab_sea12
    runsroot=/scratch/08382/shoshi/labsea_runs
    mitgcmutils=/work2/08382/shoshi/stampede3/MITgcm_c69j/MITgcm/utils/python/MITgcmutils
    optfile=/work2/08382/shoshi/stampede3/computing/optfiles/linux_amd64_ifort+mpi_stampede3_skx
    scheduler=slurm
    mpi_launch="ibrun -n"
    jra3q_src=( '/scratch/08382/shoshi/jra3q/jra3q_*_2024' )
    ;;

  pfe)
    rootdir=/nobackup/sreich/MITgcm_c69j/lab_sea12
    runsroot=/nobackup/sreich/labsea_runs
    mitgcmutils=/nobackup/sreich/MITgcm_c69j/MITgcm/utils/python/MITgcmutils
    optfile=${rootdir}/linux_amd64_ifort+mpi_ice_nas
    scheduler=pbs
    mpi_launch="mpiexec -np"
    # gpcp rain comes from our own copy; everything else from atnguye4's.
    # Order matters: ln -s will not clobber the gpcp link made first.
    jra3q_src=( '/nobackup/sreich/jra3q/jra3q_gpcp_rain_2024'
                '/nobackupp17/atnguye4/jra3q/jra3q_*_2024' )
    ;;

  *)
    echo "machine_config.sh: ERROR -- unknown LABSEA_MACHINE '${LABSEA_MACHINE}'" >&2
    echo "  Expected one of: sverdrup stampede3 pfe" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac

export rootdir runsroot mitgcmutils optfile scheduler mpi_launch

# ------------------------------------------- 3. tiling / rank-count defaults --
# These follow SIZE.h + data.exch2 of the build you point at, NOT the machine.
# nprocs = (number of tiles) - (length of blankList).
#   HR 32x32 tiling: 15x12 = 180 tiles - 77 blanks = 103   <- current default
#   LR  24x24 tiling:  5x4 =  20 tiles -  4 blanks =  16
# A driver may override these after sourcing, but 103/16 is what the committed
# SIZE.h / data.exch2 pair describes.
nprocs_hr=${nprocs_hr:-103}
nprocs_lr=${nprocs_lr:-16}

# ------------------------------------------------------------- 4. functions --

# labsea_run <nprocs> <executable> [args...]
labsea_run() {
    local n=$1; shift
    ${mpi_launch} "${n}" "$@"
}

labsea_load_modules() {
    case "${LABSEA_MACHINE}" in
      sverdrup)
        module load intel impi netcdf netcdf-fortran phdf5
        ulimit -s unlimited
        ;;
      stampede3)
        module purge
        module load intel/25.1 impi/21.15 netcdf/4.9.2 hdf5/1.14.6
        export I_MPI_DEBUG=4
        ulimit -s unlimited
        ulimit -v unlimited
        ;;
      pfe)
        module use -a /swbuild/analytix/tools/modulefiles
        module load miniconda3/v4
        export FORT_BUFFERED=1
        export MPI_IB_RAILS=2
        export MPI_BUFS_PER_PROC=128
        export MPI_DISPLAY_SETTINGS
        export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:${HOME}/lib
        umask 022
        ulimit -s hard
        ulimit -u hard
        ;;
    esac
}

# labsea_conda_activate [envname]   -- default py38
labsea_conda_activate() {
    local env=${1:-py38}
    case "${LABSEA_MACHINE}" in
      pfe)
        # miniconda3/v4 module provides `source activate`, not `conda activate`
        source activate "${env}"
        ;;
      *)
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "${env}"
        ;;
    esac
}

labsea_conda_deactivate() {
    case "${LABSEA_MACHINE}" in
      pfe)     source deactivate 2>/dev/null || conda deactivate 2>/dev/null || true ;;
      *)       conda deactivate ;;
    esac
}

# labsea_link_jra3q [destdir]   -- default ./jra3q
labsea_link_jra3q() {
    local dest=${1:-./jra3q}
    mkdir -p "${dest}"
    local pattern
    for pattern in "${jra3q_src[@]}"; do
        # unquoted on purpose: the glob must expand here
        ln -s ${pattern} "${dest}/" 2>/dev/null
    done
    return 0
}
