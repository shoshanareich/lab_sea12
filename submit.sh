#!/bin/bash
# ---------------------------------------------------------------------------
# submit.sh -- submit a driver script to whichever scheduler this machine uses.
#
#     ./submit.sh run_optimization.sh
#     ./submit.sh run_swot_optimization.sh -t 24:00:00     # extra args passed through
#     ./submit.sh -n run_optimization.sh                   # dry run: print, don't submit
#
# Detects sverdrup / stampede3 / pfe via machine_config.sh, picks the matching
# thin wrapper (submit_slurm_*.sh / submit_pbs_*.sh), and hands it the driver
# name in $DRIVER. The wrapper carries the scheduler directives; the driver
# carries the science.
#
# Override the detected machine with:
#     export LABSEA_MACHINE=sverdrup|stampede3|pfe
# ---------------------------------------------------------------------------
set -e

here=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
source "${here}/machine_config.sh"

dryrun=false
if [ "$1" = "-n" ] || [ "$1" = "--dry-run" ]; then
    dryrun=true
    shift
fi

driver=$1
if [ -z "${driver}" ]; then
    echo "usage: ./submit.sh [-n] <driver.sh> [extra scheduler args...]" >&2
    echo >&2
    echo "available drivers:" >&2
    ls "${here}"/run_*.sh 2>/dev/null | xargs -n1 basename | sed 's/^/    /' >&2
    exit 1
fi
shift

driver=$(basename "${driver}")
if [ ! -f "${here}/${driver}" ]; then
    echo "submit.sh: no such driver '${driver}' in ${here}" >&2
    exit 1
fi

job=${driver%.sh}

case "${scheduler}" in

  slurm)
    wrapper="${here}/submit_slurm_${LABSEA_MACHINE}.sh"
    [ -f "${wrapper}" ] || { echo "submit.sh: missing ${wrapper}" >&2; exit 1; }
    # Command-line args OVERRIDE in-file directives, so only supply a name /
    # output path when the wrapper does not declare one itself. Whatever you
    # put in the wrapper wins.
    naming=()
    grep -qE '^#SBATCH[[:space:]]+(-J|--job-name)' "${wrapper}" || naming+=(-J "${job}")
    grep -qE '^#SBATCH[[:space:]]+(-o|--output)'   "${wrapper}" || naming+=(-o "${here}/${job}.%j.out")
    grep -qE '^#SBATCH[[:space:]]+(-e|--error)'    "${wrapper}" || naming+=(-e "${here}/${job}.%j.err")
    set -- sbatch "${naming[@]}" \
                  --export=ALL,DRIVER="${driver}" \
                  "$@" "${wrapper}"
    ;;

  pbs)
    wrapper="${here}/submit_pbs_${LABSEA_MACHINE}.sh"
    [ -f "${wrapper}" ] || { echo "submit.sh: missing ${wrapper}" >&2; exit 1; }
    # Same rule as SLURM: the wrapper's own -N/-o/-e win if present.
    # PBS job names are length-limited, hence the truncation.
    naming=()
    grep -qE '^#PBS[[:space:]]+-N' "${wrapper}" || naming+=(-N "${job:0:15}")
    grep -qE '^#PBS[[:space:]]+-o' "${wrapper}" || naming+=(-o "${here}/${job}.out")
    grep -qE '^#PBS[[:space:]]+-e' "${wrapper}" || naming+=(-e "${here}/${job}.err")
    set -- qsub "${naming[@]}" \
                -v DRIVER="${driver}" \
                "$@" "${wrapper}"
    ;;

  *)
    echo "submit.sh: unknown scheduler '${scheduler}' for ${LABSEA_MACHINE}" >&2
    exit 1
    ;;
esac

echo "machine : ${LABSEA_MACHINE} (${scheduler})"
echo "driver  : ${driver}"
echo "command : $*"

if ${dryrun}; then
    echo "(dry run -- not submitted)"
else
    "$@"
fi
