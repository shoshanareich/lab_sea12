#!/bin/bash
#SBATCH -J labsea12 
#SBATCH -o labsea12.%j.out
#SBATCH -e labsea12.%j.err
#SBATCH -N 4
#SBATCH -n 103
#SBATCH -t 48:00:00
#SBATCH --mail-user=sreich@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
# ---------------------------------------------------------------------------
# SLURM wrapper for sverdrup.  Holds the scheduler directives only -- all the
# actual work lives in the driver script named by $DRIVER.
#
# Do not run this directly; use:
#     ./submit.sh run_optimization.sh
# which sets -J/-o/-e and exports DRIVER for you.
#
# Sizing: 28 cores/node here, so 103 HR ranks -> 4 nodes.
# ---------------------------------------------------------------------------

# NOTE: SLURM and PBS copy the batch script into a spool directory before
# running it, so BASH_SOURCE does NOT point at the repo at run time.
# Use the path submit.sh records in LABSEA_ROOT, falling back to the
# scheduler's own submit-directory variable for a hand-rolled submission.
repo="${LABSEA_ROOT:-${SLURM_SUBMIT_DIR:-}}"
if [ -z "${repo}" ] || [ ! -f "${repo}/machine_config.sh" ]; then
    echo "submit_slurm_sverdrup.sh: cannot locate the lab_sea12 checkout." >&2
    echo "  Expected machine_config.sh in: ${repo:-<unset>}" >&2
    echo "  Use ./submit.sh <driver.sh>, or export LABSEA_ROOT=<checkout>." >&2
    exit 1
fi
source "${repo}/machine_config.sh"

if [ -z "${DRIVER}" ]; then
    echo "submit_slurm_sverdrup.sh: DRIVER is not set." >&2
    echo "  Either:  ./submit.sh <driver.sh>" >&2
    echo "  or:      sbatch --export=ALL,DRIVER=run_optimization.sh submit_slurm_sverdrup.sh" >&2
    exit 1
fi

cd "${repo}" || exit 1
exec bash "./${DRIVER}"
