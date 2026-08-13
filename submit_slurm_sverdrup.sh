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

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/machine_config.sh"

if [ -z "${DRIVER}" ]; then
    echo "submit_slurm_sverdrup.sh: DRIVER is not set." >&2
    echo "  Either:  ./submit.sh <driver.sh>" >&2
    echo "  or:      sbatch --export=ALL,DRIVER=run_optimization.sh submit_slurm_sverdrup.sh" >&2
    exit 1
fi

cd "${rootdir}" || exit 1
exec bash "./${DRIVER}"
