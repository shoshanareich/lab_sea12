#!/bin/bash
#PBS -q long
#PBS -l select=3:ncpus=40:model=sky_ele
#PBS -l walltime=48:00:00
#PBS -N argoMG
#PBS -o argoMG.out
#PBS -e argoMG.err
#PBS -W umask=33
#PBS -m bea
# ---------------------------------------------------------------------------
# PBS wrapper for NASA Pleiades.  Holds the scheduler directives only -- all
# the actual work lives in the driver script named by $DRIVER.
#
# Do not run this directly; use:
#     ./submit.sh run_optimization.sh
# which sets -N/-o/-e and passes DRIVER through -v for you.
#
# Sizing: sky_ele is 40 cores/node, so 103 HR ranks -> 3 nodes.
# ---------------------------------------------------------------------------

# NOTE: SLURM and PBS copy the batch script into a spool directory before
# running it, so BASH_SOURCE does NOT point at the repo at run time.
# Use the path submit.sh records in LABSEA_ROOT, falling back to the
# scheduler's own submit-directory variable for a hand-rolled submission.
repo="${LABSEA_ROOT:-${PBS_O_WORKDIR:-}}"
if [ -z "${repo}" ] || [ ! -f "${repo}/machine_config.sh" ]; then
    echo "submit_pbs_pfe.sh: cannot locate the lab_sea12 checkout." >&2
    echo "  Expected machine_config.sh in: ${repo:-<unset>}" >&2
    echo "  Use ./submit.sh <driver.sh>, or export LABSEA_ROOT=<checkout>." >&2
    exit 1
fi
source "${repo}/machine_config.sh"

if [ -z "${DRIVER}" ]; then
    echo "submit_pbs_pfe.sh: DRIVER is not set." >&2
    echo "  Either:  ./submit.sh <driver.sh>" >&2
    echo "  or:      qsub -v DRIVER=run_optimization.sh submit_pbs_pfe.sh" >&2
    exit 1
fi

cd "${repo}" || exit 1
exec bash "./${DRIVER}"
