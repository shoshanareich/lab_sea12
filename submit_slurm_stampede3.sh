#!/bin/bash
#SBATCH -J ArgoMG
#SBATCH -o ArgoMG.%j.out
#SBATCH -e ArgoMG.%j.err
#SBATCH -t 48:00:00
#SBATCH -p skx
#SBATCH -N 6
#SBATCH -n 103
#SBATCH -A OCE23001
#SBATCH --mail-user=sreich@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
# ---------------------------------------------------------------------------
# SLURM wrapper for TACC Stampede3.  Holds the scheduler directives only --
# all the actual work lives in the driver script named by $DRIVER.
#
# Do not run this directly; use:
#     ./submit.sh run_optimization.sh
# which sets -J/-o/-e and exports DRIVER for you.
# ---------------------------------------------------------------------------

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/machine_config.sh"

if [ -z "${DRIVER}" ]; then
    echo "submit_slurm_stampede3.sh: DRIVER is not set." >&2
    echo "  Either:  ./submit.sh <driver.sh>" >&2
    echo "  or:      sbatch --export=ALL,DRIVER=run_optimization.sh submit_slurm_stampede3.sh" >&2
    exit 1
fi

cd "${rootdir}" || exit 1
exec bash "./${DRIVER}"
