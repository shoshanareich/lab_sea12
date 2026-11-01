#!/usr/bin/env python3
"""
Poster version of mg_27: along-track SSH spectra (top) and coherence^2 with
SWOT (bottom), for the 10-day sigmaL runs at iteration 10 in both.

    SWOT observations
    model iter 0           shared by both runs
    multi-grid  it0010     assim_swot_MG_10days_sigmaL (x1 LR viscosity), accepted
    HR-only     it0010     assim_swot_HR_10days_sigmaL, accepted, its last iterate

Differences from plot_swot_spectra_HR_vs_MG.py --variant coherence:
  * obs file is the _sigmaL one.  Its obs_val is identical to the old file's
    except for 3 samples set to -9999 (excluded from the cost, section 17 of
    the debug notes); those are masked here -- the old code tests only isfinite.
  * one legend, in the top panel, covering both panels.
  * poster title and font sizes; PDF output in figs/.

All method notes (coherence as a ratio of sums, 15 log bins, MIN_NEFF bins
dropped, SSH not SSHA so lambda > 200 km is MDT) are in the parent script.

Usage:  python3 plot_swot_coherence_poster_sigmaL.py
"""

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import plot_swot_spectra as pss
import plot_swot_spectra_HR_vs_MG as hm

RUNSROOT = pss.RUNSROOT
OBS_NAME = 'swot_obsfit_cycles_9thru11_labsea_L3v3_sigmaL'
OBS_FILE = '/home/shoshi/MITgcm_c69j/lab_sea12/input_binaries_hires/' + OBS_NAME + '.nc'
OUTFILE  = '/home/shoshi/MITgcm_c69j/lab_sea12/figs/swot_spectra_coherence_HR_vs_MG_sigmaL_iter10.pdf'

HR = RUNSROOT + 'assim_swot_HR_10days_sigmaL/'
MG = RUNSROOT + 'assim_swot_MG_10days_sigmaL/'

# (label, run dir, iteration, colour, linewidth, linestyle)
CURVES = [
    ('Model iter 0',        MG,  0, hm.CRIMSON, 1.8, ':'),
    ('Multi-grid  iter 10', MG, 10, hm.BLUE,    2.4, '-'),
    ('HR-only  iter 10',    HR, 10, hm.GREEN,   2.4, '-'),
]

TITLE = 'Along-track SSH in the Labrador Sea: SWOT vs. assimilation, 10-day window'


def main():
    plt.rcParams.update({'font.size': 14, 'pdf.fonttype': 42})

    ds = xr.open_dataset(OBS_FILE)
    nobs = ds.dims['iOBS']
    obs = ds.obs_val.values.astype(float)
    obs[obs <= -9999.0] = np.nan

    loaded = []
    for label, run, it, color, lw, ls in CURVES:
        mv, mm = pss.load_equi(run + 'run_adhi_it%04d/' % it, nobs, obs_name=OBS_NAME)
        loaded.append((label, mv, mm, color, lw, ls))

    usable = np.isfinite(obs)
    for _, mv, mm, _, _, _ in loaded:
        usable &= (mm == 1) & np.isfinite(mv)
    segments = pss.build_segments(ds, usable)
    r_swath, n_swaths = hm.swath_redundancy(ds, segments)
    print('usable samples: %d of %d' % (usable.sum(), nobs))
    print('segments: %d  (%d points, %d swaths)'
          % (len(segments), sum(len(s) for s in segments), n_swaths))

    # top: PSD on the 120-bin interp grid, as mg_27
    pss.set_grid(120, 'interp')
    f_psd = pss.master_frequencies
    mean_o, ci_o = pss.stack(segments, obs)[:2]
    psd = [(label, pss.stack(segments, mv), color, lw, ls)
           for label, mv, mm, color, lw, ls in loaded]

    # bottom: coherence on 15 log bins, bins with n_eff < MIN_NEFF dropped
    nb = hm.N_BINS_COH
    pss.set_grid(nb, 'logbin')
    f_coh, edges = pss.master_frequencies, pss._edges
    coh, n_d = [], None
    for label, mv, mm, color, lw, ls in loaded:
        c2, nd = hm.coherence(segments, obs, mv, edges, nb)
        n_d = nd if n_d is None else n_d
        coh.append((label, c2, color, lw, ls))
    keep = n_d * r_swath >= hm.MIN_NEFF
    coh = [(lb, np.where(keep, c2, np.nan), c, lw, ls) for lb, c2, c, lw, ls in coh]

    lams = (1000, 500, 200, 150, 100, 70, 50, 35)
    print('\ncoh^2   ' + ''.join('%-8d' % l for l in lams))
    for label, c2, _, _, _ in coh:
        i = [np.argmin(np.abs(1 / f_coh - l)) for l in lams]
        print('%-22s' % label + ''.join('%-8.3f' % c2[j] for j in i))

    # ------------------------------------------------------------------ plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10.4), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1.55, 1],
                                                    hspace=0.08))

    ax1.loglog(f_psd, mean_o, color='black', lw=2.6, label='SWOT observations')
    ax1.fill_between(f_psd, np.maximum(mean_o - ci_o, 1e-10), mean_o + ci_o,
                     color='black', alpha=0.15)
    for label, (mean, ci, _, _), color, lw, ls in psd:
        ax1.loglog(f_psd, mean, color=color, lw=lw, ls=ls, label=label)
        ax1.fill_between(f_psd, np.maximum(mean - ci, 1e-10), mean + ci,
                         color=color, alpha=0.15)
    for ax in (ax1, ax2):
        ax.axvspan(1 / hm.MODEL_RES, pss.F_NYQ, color='0.85', alpha=0.5, zorder=0)
        ax.grid(True, which='both', ls='--', alpha=0.35)
    ax1.text(0.033, 0.05, 'below model\neffective res.', fontsize=11,
             color='0.3', ha='center', transform=ax1.get_xaxis_transform())
    ax1.set_ylabel('PSD (m$^2$ / cycles km$^{-1}$)')
    ax1.set_ylim(1e-6, 2e2)
    ax1.legend(loc='lower left', fontsize=13, framealpha=0.95)
    ax1.set_title(TITLE, pad=16, fontsize=15, fontweight='bold')

    axt = ax1.secondary_xaxis('top', functions=(lambda x: 1 / x, lambda x: 1 / x))
    axt.set_xlabel('Wavelength (km)')
    axt.set_xticks([2000, 1000, 500, 200, 100, 50, 30, 20])
    axt.set_xticklabels(['2000', '1000', '500', '200', '100', '50', '30', '20'])

    for label, c2, color, lw, ls in coh:
        ax2.semilogx(f_coh, c2, color=color, lw=lw, ls=ls, marker='o', ms=5, zorder=3)
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel('Coherence$^2$ with SWOT')
    ax2.set_xlabel('Wavenumber (cycles km$^{-1}$)')

    fig.tight_layout()
    fig.savefig(OUTFILE)
    print('\nwrote %s' % OUTFILE)


if __name__ == '__main__':
    main()
