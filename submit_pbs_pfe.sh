#!/bin/bash
#PBS -q long
#PBS -l select=3:ncpus=40:model=sky_ele
#PBS -l walltime=48:00:00
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

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/machine_config.sh"

if [ -z "${DRIVER}" ]; then
    echo "submit_pbs_pfe.sh: DRIVER is not set." >&2
    echo "  Either:  ./submit.sh <driver.sh>" >&2
    echo "  or:      qsub -v DRIVER=run_optimization.sh submit_pbs_pfe.sh" >&2
    exit 1
fi

cd "${rootdir}" || exit 1
exec bash "./${DRIVER}"
