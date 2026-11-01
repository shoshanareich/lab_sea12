#!/usr/bin/env python3
"""
Increment spectrum: PSD of (mod_it - mod_it0), the spectrum of what the
assimilation actually changed in SSH along the SWOT tracks.

This is the direct test of "does the coarse adjoint cap the scale of the
correction?".  MG and HR share a 1/12 deg forward model and are bit-identical
at iteration 0, so they differ in exactly one thing: adjoint resolution
(1/3 deg vs 1/12 deg).  If that resolution capped the correction scale, MG
would have to fall off FASTER at high wavenumber than HR -- systematically, and
under both control priors.  LR is shown for contrast: it changes the forward
model too, and its collapse below 200 km is a forward-model limit, not an
adjoint one.

Left panel normalises each curve at 400 km, so it compares the SHAPE of the
correction independent of how large a step the optimizer took.  Right panel is
absolute.

Iterations default to the last ACCEPTED m1qn3 iterate of each run -- m1qn3
evaluates the cost at trial points too, so it#### is not automatically an
iterate (MG v1 it5, used by plot_misfit_scale.py's panel (d), is a rejected
trial).  Verify with plot_mg_diagnostics.classify_iterates before changing them.

Usage:  python3 plot_increment_spectra.py [outfile.png]
"""

import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/shoshi/MITgcm_c69j/lab_sea12/mappings')
from plot_swot_spectra import (load_equi, build_segments, stack, set_grid,
                               OBS_FILE, RUNSROOT)
import plot_swot_spectra as P

# ============================== configuration ==============================
# label: (run dir, iteration, colour, linestyle, linewidth)
CURVES = [
    ('MG v1  iter 8',   'assim_swot_MG_10days_xxmult16_newCtrlW/',     8, 'royalblue', '--', 1.7),
    ('MG v2  iter 5',   'assim_swot_MG_10days_xxmult16_newCtrlW_v2/',  5, 'royalblue', '-',  2.0),
    ('HR v1  iter 17',  'assim_swot_HR_10days_newCtrlW/',             17, '#1baf7a',   '--', 1.7),
    ('HR v2  iter 10',  'assim_swot_HR_10days_newCtrlW_v2/',          10, '#1baf7a',   '-',  2.0),
    ('LR only  iter 5', 'assim_swot_LR_10days/',                       5, '#eda100',   '-.', 1.7),
]
SUBPAT   = {'assim_swot_LR_10days/': 'run_adlo_it%04d/'}    # LR writes run_adlo
N_BINS   = 40
NORM_LAM = 400.0
LAM_MAX  = 700.0
MODEL_RES = 60.0
# ===========================================================================

set_grid(N_BINS, 'logbin')
F = P.master_frequencies


def at(S, lam):
    return S[np.argmin(np.abs(1 / F - lam))]


def main(outfile='increment_spectra.png'):
    ds = xr.open_dataset(OBS_FILE)
    nobs = ds.dims['iOBS']

    mod, usable = {}, np.isfinite(ds.obs_val.values)
    for label, run, it, _, _, _ in CURVES:
        pat = SUBPAT.get(run, 'run_adhi_it%04d/')
        for k in (0, it):
            mv, mm = load_equi(RUNSROOT + run + pat % k, nobs)
            mod[(run, k)] = mv
            usable &= (mm == 1) & np.isfinite(mv)
    segments = build_segments(ds, usable)
    print('%d usable samples, %d segments' % (usable.sum(), len(segments)))

    inc = {label: stack(segments, mod[(run, it)] - mod[(run, 0)])[0]
           for label, run, it, _, _, _ in CURVES}

    lams = [400, 300, 200, 150, 100, 70, 50, 35, 25]
    for title, norm in (('normalised at %g km' % NORM_LAM, True), ('absolute', False)):
        print('\n=== increment PSD, %s ===' % title)
        print('%-16s ' % 'lam(km)' + ''.join('%-10d' % l for l in lams))
        for label, *_ in CURVES:
            S = inc[label]
            v = [at(S, l) / (at(S, NORM_LAM) if norm else 1.0) for l in lams]
            fmt = '%-10.3f' if norm else '%-10.2e'
            print('%-16s ' % label + ''.join(fmt % x for x in v))

    # MG-vs-HR difference, within each prior
    print('\n=== MG / HR increment ratio, within each prior (normalised shape) ===')
    print('%-16s ' % 'lam(km)' + ''.join('%-10d' % l for l in lams))
    for ver, m, h in (('v1', 'MG v1  iter 8', 'HR v1  iter 17'),
                      ('v2', 'MG v2  iter 5', 'HR v2  iter 10')):
        r = [(at(inc[m], l) / at(inc[m], NORM_LAM)) /
             (at(inc[h], l) / at(inc[h], NORM_LAM)) for l in lams]
        print('%-16s ' % ('MG/HR ' + ver) + ''.join('%-10.2f' % x for x in r))

    fig, ax = plt.subplots(1, 2, figsize=(14.5, 6.2))
    for j, norm in enumerate((True, False)):
        a = ax[j]
        for label, run, it, color, ls, lw in CURVES:
            S = inc[label]
            y = S / at(S, NORM_LAM) if norm else S
            a.loglog(F, y, color=color, ls=ls, lw=lw, marker='o', ms=3, label=label)
        a.axvspan(1 / MODEL_RES, 1 / (2 * 10.0), color='0.85', alpha=0.5, zorder=0)
        a.set_xlim(1 / LAM_MAX, 5.5e-2)
        a.grid(True, which='both', ls='--', alpha=0.35)
        a.set_xlabel('Wavenumber (cycles / km)')
        a.set_ylabel('increment PSD, normalised at %g km' % NORM_LAM if norm
                     else 'increment PSD (m$^2$ / cycles/km)')
        a.set_title('(%s) %s' % ('ab'[j], 'shape of the correction'
                    if norm else 'absolute size of the correction'),
                    loc='left', fontsize=11)
        t = a.secondary_xaxis('top', functions=(lambda x: 1/x, lambda x: 1/x))
        t.set_xticks([500, 200, 100, 50, 25])
        t.set_xticklabels(['500', '200', '100', '50', '25'])
        t.set_xlabel('Wavelength (km)')
        if norm:
            a.set_ylim(1e-4, 2); a.legend(fontsize=8.5, loc='lower left')

    fig.suptitle('What the assimilation actually changed:  PSD(mod$_{it}$ - mod$_{it0}$)\n'
                 'MG and HR share a 1/12$^\\circ$ forward model, so they differ only in '
                 'adjoint resolution', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(outfile, dpi=155)
    print('\nwrote %s' % outfile)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'increment_spectra.png')
