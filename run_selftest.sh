#!/bin/bash
# ---------------------------------------------------------------------------
# run_selftest.sh
#
# A throwaway driver that exercises the whole submission chain in ~1 minute
# without running the model:
#
#   submit.sh -> scheduler -> spooled wrapper -> LABSEA_ROOT resolution
#             -> machine_config.sh -> labsea_load_modules -> labsea_run (MPI)
#             -> labsea_conda_activate -> mappings/machine_config.py
#
# This is the part that dry runs CANNOT reach: the wrapper only executes for
# real when the scheduler runs it, from its spool directory.
#
#     ./submit.sh run_selftest.sh -t 00:05:00
#
# Reports PASS/FAIL per check and exits non-zero if anything failed.
# ---------------------------------------------------------------------------

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/machine_config.sh"

fails=0
ok()   { printf '  PASS  %s\n' "$*"; }
bad()  { printf '  FAIL  %s\n' "$*"; fails=$(( fails + 1 )); }
note() { printf '  ..    %s\n' "$*"; }

echo "==================================================================="
echo " lab_sea12 submission self-test"
echo " host          : $(hostname -f)"
echo " machine       : ${LABSEA_MACHINE}"
echo " scheduler     : ${scheduler}"
echo " cwd           : $(pwd)"
echo " LABSEA_ROOT   : ${LABSEA_ROOT:-<unset - not launched via submit.sh>}"
echo " SLURM job     : ${SLURM_JOB_ID:-none}"
echo " PBS job       : ${PBS_JOBID:-none}"
echo "==================================================================="

# --- 1. did the wrapper hand us a sane repo? -------------------------------
echo
echo "[1] repo resolution"
[ -f "${rootdir}/machine_config.sh" ] && ok "machine_config.sh found under rootdir" \
                                      || bad "machine_config.sh NOT under rootdir=${rootdir}"
# Compare resolved paths: on pfe /nobackup/sreich is a symlink to
# /nobackupp18/sreich, so cwd and rootdir are the same place under two names.
if [ "$(readlink -f "$(pwd)")" = "$(readlink -f "${rootdir}")" ]; then
    ok "cwd resolves to rootdir"
    [ "$(pwd)" != "${rootdir}" ] && note "(via symlink: $(pwd) -> ${rootdir})"
else
    note "cwd $(pwd) != rootdir ${rootdir} (fine if you ran this by hand)"
fi
[ -d "${runsroot}" ] && ok "runsroot exists: ${runsroot}" \
                     || bad "runsroot missing: ${runsroot}"

# --- 2. modules ------------------------------------------------------------
# IMPORTANT: do NOT pipe labsea_load_modules -- a pipeline runs it in a
# subshell and every `module load` is discarded, so later checks would run
# without the environment they need. Redirect to a file instead.
echo
echo "[2] labsea_load_modules"
modlog=$(mktemp "${TMPDIR:-/tmp}/labsea_mod.XXXXXX")
if labsea_load_modules > "${modlog}" 2>&1; then
    ok "labsea_load_modules returned cleanly"
else
    bad "labsea_load_modules failed"
fi
sed 's/^/        /' "${modlog}" | head -20
rm -f "${modlog}"
note "loaded: $( { module list; } 2>&1 | tr '\n' ' ' | sed 's/  */ /g' | cut -c1-240 )"

# --- 3. the executables the drivers expect ---------------------------------
echo
echo "[3] builds"
for b in build_adhi_2lev_seaice_update_mpi build_adlo_2lev_seaice_update_mpi; do
    if [ -x "${rootdir}/${b}/mitgcmuv_ad" ]; then
        ok "${b}/mitgcmuv_ad"
    else
        bad "${b}/mitgcmuv_ad missing or not executable"
    fi
done

# --- 4. the MPI launcher ---------------------------------------------------
# The single most valuable check here: ibrun vs mpiexec, and whether the
# allocation actually gives us ranks.
echo
echo "[4] labsea_run (mpi_launch = '${mpi_launch}')"
note "launcher: $(command -v "${mpi_launch%% *}" 2>/dev/null || echo '<NOT ON PATH>')"
if out=$(labsea_run 2 /bin/hostname 2>&1); then
    ok "launched 2 ranks"
    echo "${out}" | sed 's/^/        /' | head -4
else
    bad "labsea_run failed:"
    echo "${out}" | sed 's/^/        /' | head -6
    # HPE MPT (pfe) can refuse to launch a binary that never calls MPI_Init,
    # so a failure here on a serial command is not proof that the real model
    # launch is broken. Retry with the LR executable's --help-ish no-op is not
    # safe, so just report enough to tell the two cases apart.
    note "if this says 'could not run executable' on pfe, check whether an MPI"
    note "module (e.g. mpi-hpe/mpt) is loaded -- labsea_load_modules does not"
    note "load one, and confirm with a real mitgcmuv_ad launch before trusting it."
fi

# --- 5. conda + the Python side of the config ------------------------------
echo
echo "[5] labsea_conda_activate + mappings/machine_config.py"
if labsea_conda_activate 2>/dev/null; then
    ok "conda env activated ($(python3 --version 2>&1))"
    if out=$(cd "${rootdir}/mappings" && python3 -c "
import machine_config as m
m.add_mitgcmutils()
from MITgcmutils import rdmds
print('machine=%s rootdir=%s' % (m.MACHINE, m.ROOTDIR))
print('MITgcmutils OK')
" 2>&1); then
        ok "machine_config.py + MITgcmutils import"
        echo "${out}" | sed 's/^/        /'
    else
        bad "python import failed:"
        echo "${out}" | tail -4 | sed 's/^/        /'
    fi
    labsea_conda_deactivate 2>/dev/null
else
    bad "labsea_conda_activate failed"
    note "conda   : $(command -v conda 2>/dev/null || echo '<not on PATH>')"
    note "activate: $(command -v activate 2>/dev/null || echo '<not on PATH>')"
    note "on pfe this comes from 'module load miniconda3/v4' in labsea_load_modules"
fi

# --- 6. grid dirs the mapping scripts read ---------------------------------
echo
echo "[6] grid directories under runsroot"
for g in grid_hires grid_lores grid_lores_cleanbathy grid_lores_cleanbathy_v2 grid_hires_cleanbathy; do
    [ -d "${runsroot}/${g}" ] && ok "${g}" || bad "${g} missing under ${runsroot}"
done

# --- verdict ---------------------------------------------------------------
echo
echo "==================================================================="
if [ ${fails} -eq 0 ]; then
    echo " SELF-TEST PASSED on ${LABSEA_MACHINE}"
else
    echo " SELF-TEST: ${fails} FAILURE(S) on ${LABSEA_MACHINE}"
fi
echo "==================================================================="
exit $(( fails > 0 ))
