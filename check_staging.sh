#!/bin/bash
# ---------------------------------------------------------------------------
# check_staging.sh -- dry-run the link_*.sh scripts without touching anything.
#
# Runs each link script for real, in a throwaway sandbox, with cp / ln / rm /
# sed / mv shadowed by logging stubs. So it exercises the actual control flow
# (the iter<1 vs iter>0 branches, labsea_link_jra3q, the namelist seds) and
# reports every file each script WOULD stage -- and which of those are missing.
#
#   ./check_staging.sh                  # check this machine (auto-detected)
#   ./check_staging.sh -m pfe           # form paths as pfe would (see note below)
#   ./check_staging.sh -v               # also list the sources that were found
#   ./check_staging.sh -s link_lores.sh # just one script
#
# NOTE ON -m: you can only stat files that exist on the machine you are sitting
# on. Run with -m for another machine and it checks path FORMATION only -- that
# every variable expanded, nothing came out empty or with a literal ${...} left
# in it. To truly verify another machine's staging, copy this script there (it
# needs no scheduler and no allocation) and run it bare on the login node.
# ---------------------------------------------------------------------------
set -u

here=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")

verbose=false
only=""
while [ $# -gt 0 ]; do
    case "$1" in
      -m|--machine) export LABSEA_MACHINE="$2"; shift 2 ;;
      -v|--verbose) verbose=true; shift ;;
      -s|--script)  only="$2"; shift 2 ;;
      -h|--help)    sed -n '2,25p' "${BASH_SOURCE[0]}"; exit 0 ;;
      *) echo "check_staging.sh: unknown option '$1'" >&2; exit 1 ;;
    esac
done

export LABSEA_QUIET=1
source "${here}/machine_config.sh"

# Is the config we are testing the machine we are actually sitting on?
native=false
[ -d "${rootdir}" ] && native=true

echo "=============================================================="
echo " machine   : ${LABSEA_MACHINE}"
echo " rootdir   : ${rootdir}"
echo " runsroot  : ${runsroot}"
if ${native}; then
    echo " mode      : FULL -- rootdir exists here, so files will be stat'ed"
else
    echo " mode      : FORMATION ONLY -- rootdir does not exist on this host."
    echo "             Checking that paths expand, not that files are present."
fi
echo "=============================================================="

sandbox=$(mktemp -d "${TMPDIR:-/tmp}/labsea_staging.XXXXXX")
trap 'rm -rf "${sandbox}"' EXIT

log="${sandbox}/staged.log"
: > "${log}"

# --- the stubs -------------------------------------------------------------
# Record sources instead of copying/linking. Everything else is a no-op so the
# script runs to completion.
_record() {
    local kind=$1; shift
    # drop trailing destination arg; log each remaining source
    local args=("$@") n=$(( $# - 1 )) i
    for (( i=0; i<n; i++ )); do
        case "${args[i]}" in
          -*) continue ;;                       # flags
        esac
        printf '%s\t%s\n' "${kind}" "${args[i]}" >> "${log}"
    done
}
cp()  { _record cp "$@"; }
ln()  { _record ln "$@"; }
rm()  { :; }
sed() { :; }
mv()  { :; }
export -f _record cp ln rm sed 2>/dev/null || true

run_one() {
    local script=$1; shift
    echo
    echo "### ${script}  $*"
    : > "${log}"
    (
        cd "${sandbox}" || exit 1
        # shellcheck disable=SC1090
        source "${here}/${script}" "$@"
    ) >/dev/null 2>&1

    local total=0 miss=0 prod=0
    while IFS=$'\t' read -r kind src; do
        [ -z "${src}" ] && continue
        total=$(( total + 1 ))
        # classify: repo input (under rootdir) vs run product (under runsroot)
        local class="input"
        case "${src}" in
          "${runsroot}"*) class="product" ;;
        esac
        # unexpanded or empty variable is always an error, on any machine
        if [ -z "${src}" ] || [[ "${src}" == *'${'* ]] || [[ "${src}" == //* ]]; then
            echo "  BAD-EXPANSION  ${src}"; miss=$(( miss + 1 )); continue
        fi
        if ${native} && [ "${class}" = "input" ]; then
            if [ "$(eval "ls -d ${src}" 2>/dev/null | wc -l)" -eq 0 ]; then
                echo "  MISSING        ${src}"; miss=$(( miss + 1 ))
            elif ${verbose}; then
                echo "  ok             ${src}"
            fi
        else
            [ "${class}" = "product" ] && prod=$(( prod + 1 ))
            ${verbose} && echo "  ${class}        ${src}"
        fi
    done < "${log}"

    echo "  -- ${total} sources; ${miss} problem(s); ${prod} are run products (created by earlier stages)"
    return $(( miss > 0 ))
}

# --- representative arguments ---------------------------------------------
scratch="${runsroot}/assim_argo_MG"
bhi="${rootdir}/build_adhi_2lev_seaice_update_mpi"
blo="${rootdir}/build_adlo_2lev_seaice_update_mpi"

rc=0
should_run() { [ -z "${only}" ] || [ "${only}" = "$1" ]; }

# link_hires.sh <iter> <optimext> <scratchdir> <builddir>
should_run link_hires.sh && { run_one link_hires.sh 0 0000 "${scratch}" "${bhi}" || rc=1; }
should_run link_hires.sh && { run_one link_hires.sh 1 0001 "${scratch}" "${bhi}" || rc=1; }
# link_lores.sh <iter> <dirhires> <scratchdir> <builddir>
should_run link_lores.sh && { run_one link_lores.sh 1 "${scratch}/run_adhi_it0001" "${scratch}" "${blo}" || rc=1; }
# link_packunpack.sh <iter> <workdir> <optimdir> <builddir> <dirhires>
should_run link_packunpack.sh && { run_one link_packunpack.sh 1 "${scratch}/run_adlo_packunpack" \
                                   "${scratch}/OPTIM" "${blo}" "${scratch}/run_adhi_it0001" || rc=1; }

echo
if [ ${rc} -eq 0 ]; then
    echo "RESULT: no problems found for ${LABSEA_MACHINE}."
else
    echo "RESULT: problems found for ${LABSEA_MACHINE} (see MISSING / BAD-EXPANSION above)."
fi
exit ${rc}
