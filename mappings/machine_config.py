"""
machine_config.py -- shared per-machine configuration for the lab_sea12
offline mapping scripts.

Use at the top of any mappings/*.py, before importing MITgcmutils:

    from machine_config import MACHINE, ROOTDIR, RUNSROOT, add_mitgcmutils
    add_mitgcmutils()
    from MITgcmutils import rdmds

Auto-detects sverdrup / stampede3 / pfe.  To override (e.g. on an odd compute
node, or to test another machine's paths):

    export LABSEA_MACHINE=sverdrup|stampede3|pfe

This is the Python twin of ../machine_config.sh -- keep the two in sync.
"""

import os
import socket
import sys

# ------------------------------------------------------------------ detect --
def detect_machine():
    """Return ('sverdrup'|'stampede3'|'pfe', reason)."""
    override = os.environ.get('LABSEA_MACHINE')
    if override:
        if override not in CONFIG:
            raise RuntimeError(
                "LABSEA_MACHINE='%s' is not one of %s"
                % (override, ', '.join(sorted(CONFIG)))
            )
        return override, 'LABSEA_MACHINE override'

    if os.path.isdir('/nobackup/sreich'):
        return 'pfe', 'found /nobackup/sreich'
    if os.environ.get('TACC_SYSTEM') or os.path.isdir('/work2/08382/shoshi'):
        return 'stampede3', 'TACC_SYSTEM=%s / found /work2/08382/shoshi' % (
            os.environ.get('TACC_SYSTEM', 'unset'))
    if os.path.isdir('/scratch/shoshi') or socket.getfqdn().startswith('sverdrup'):
        return 'sverdrup', 'found /scratch/shoshi or hostname %s' % socket.gethostname()

    raise RuntimeError(
        "machine_config.py: cannot identify this machine. "
        "Set LABSEA_MACHINE=sverdrup|stampede3|pfe and re-run."
    )


# -------------------------------------------------------- machine settings --
CONFIG = {
    'sverdrup': dict(
        rootdir='/home/shoshi/MITgcm_c69j/lab_sea12',
        runsroot='/scratch/shoshi/labsea_MG_12',
        mitgcmutils='/home/shoshi/MITgcm_c69j/MITgcm/utils/python/MITgcmutils',
    ),
    'stampede3': dict(
        rootdir='/work2/08382/shoshi/stampede3/MITgcm_c69j/lab_sea12',
        runsroot='/scratch/08382/shoshi/labsea_runs',
        mitgcmutils='/work2/08382/shoshi/stampede3/MITgcm_c69j/MITgcm/utils/python/MITgcmutils',
    ),
    'pfe': dict(
        rootdir='/nobackup/sreich/MITgcm_c69j/lab_sea12',
        runsroot='/nobackup/sreich/labsea_runs',
        mitgcmutils='/nobackup/sreich/MITgcm_c69j/MITgcm/utils/python/MITgcmutils',
    ),
}

MACHINE, _HOW = detect_machine()

# Announce which machine we resolved to, and why -- the twin of the
# [machine_config] line printed by ../machine_config.sh.
# Silence it with: export LABSEA_QUIET=1
if not os.environ.get('LABSEA_QUIET'):
    print('[machine_config] LABSEA_MACHINE=%s  (%s)' % (MACHINE, _HOW))

ROOTDIR = CONFIG[MACHINE]['rootdir']
RUNSROOT = CONFIG[MACHINE]['runsroot']
MITGCMUTILS = CONFIG[MACHINE]['mitgcmutils']

# `dirroot` / `root_dir` in the older scripts meant the runs root, and they all
# assumed a trailing slash. Keep that alias so path concatenation still works.
DIRROOT = RUNSROOT.rstrip('/') + '/'


def add_mitgcmutils():
    """Put this machine's MITgcmutils on sys.path (idempotent)."""
    if MITGCMUTILS not in sys.path:
        sys.path.append(MITGCMUTILS)
    return MITGCMUTILS
