#!/usr/bin/env python3
"""
mg_20 -- along-track SWOT SSH wavenumber spectra, ctrlWyearly runs.
mg_27 -- the same runs' coherence with SWOT  (--variant coherence).

Same figure as mg_12b / mg_16b (interp onto a 120-bin master grid), with the
curve set the comparison needs here:

    SWOT observations
    model iter 0          identical in both runs (verified: max |mod_val diff|
                          = 0 over all 27960 samples with mod_mask = 1)
    HR-only  it0008       m1qn3 iter 8,  cost 3.923e5
    HR-only  it0020       m1qn3 iter 18, cost 2.068e5
    multi-grid it0008     cost 3.904e5

Hue carries the configuration (green = HR-only, blue = multi-grid), dashing
carries the iteration, so the two comparisons read off separate channels.

NOTE on 'MG iter 8'.  In the multi-grid run it0008 is a REJECTED line-search
trial (simul 9; m1qn3 iterates are simuls 1-6, 11-14, 18, 19, 21).  It is
plotted because it is the like-for-like same-iteration-count comparison against
HR it0008, but it is not an m1qn3 iterate.  --variant accepted swaps it for
it0011, which IS m1qn3 iter 8, and adds MG it0020 for the endpoint.

All method caveats from plot_swot_spectra.py apply unchanged -- obs_val is SSH
not SSHA so lambda > 200 km is mean dynamic topography; Nyquist is 20 km because
the obs file was smoothed and decimated; below ~60 km the model is
dissipation-limited -- though note MODEL_RES below, which is measured off the
grid rather than inherited from plot_swot_spectra's 60 km.

--variant best is the two-curve version: the lowest-cost iterate of each run,
which for both is it0020 -- HR f = 2.0676e5 (m1qn3 iter 18), MG f = 3.3034e5
(m1qn3 iter 13).  Neither run has a cost minimum in its interior, so 'best' and
'last accepted' coincide here; re-check that if either run is extended.


--variant coherence  (mg_27)
----------------------------
Two stacked panels on a shared wavenumber axis: the mg_12b PSD above, the
magnitude-squared coherence with SWOT below.  PSD says whether the model holds
the right AMOUNT of variance at a scale; coherence says whether that variance is
in the right PLACE.  Only the second one can fall while the first improves,
which is the MG failure mode recorded in section 9 of OPTIMIZATION_DEBUG_NOTES.

Three things about this panel that the spectra panel above it does not have to
worry about:

1. Coherence is a ratio of SUMS.  Sxx, Syy and Sxy must each be accumulated over
   segments and over the native frequencies in a bin BEFORE the ratio is formed;
   a per-segment, per-frequency coherence is identically 1.  That rules out the
   'interp' path the PSD panel uses, which interpolates log(PSD) per segment --
   there is no meaningful log-interpolation of a complex cross-spectrum.  The
   coherence panel therefore always log-bin averages.

2. It needs far fewer bins.  On mg_12b's 120-bin grid most bins receive one or
   two native frequencies per segment, and coh^2 from n complex products is
   biased to 1 as n -> 1 no matter what the true coherence is.  N_BINS_COH = 15
   over the same F_MIN..F_NYQ range is the coarsest grid that still reads as a
   curve against wavenumber.

3. The bias floor is drawn, because the answer is near it.  Under a true
   coherence of zero, E[coh^2] = 1/n_d for n_d independent complex products.
   Two floors are shown:

     dotted line   1/n_d       n_d = every (segment, native frequency) pair in
                               the bin, i.e. the floor IF the segments were
                               independent.  They are not.
     grey shading  1/n_eff     n_eff = n_d * n_swaths/n_segments.  track_id 0-9
                               are ten parallel lines only 10 km apart in one
                               swath, so this collapses each (cycle, pass_no)
                               swath to a single realization -- deliberately
                               pessimistic, the honest bracket being between the
                               two.

   Neither floor is DRAWN any more -- both are printed to stdout, along with
   the wavelength at which each curve falls through 1/n_eff, which is the
   number to quote.  Bins with n_eff < MIN_NEFF are still dropped rather than
   plotted at their bias.

The coherence curve set is iter 0 (shared baseline, crimson) plus MG it0005
(simul 6 = m1qn3 iter 6, f = 3.8716e5) and HR it0020 (simul 21 = m1qn3 iter 18,
f = 2.0676e5).  Note this is not best-vs-best: MG's own lowest cost is it0020 at
3.3034e5.  Colours and dashing follow the spectra figures unchanged.

Because obs_val is SSH rather than SSHA, coherence above ~200 km is coherence
with the mean dynamic topography, which any run reproduces -- read the panel
below 200 km.  (That 200 km line is not drawn.)

Usage:  python3 plot_swot_spectra_HR_vs_MG.py [outfile.png]
                 [--variant requested|accepted|best|coherence] [-n 120]
                 [-m interp|logbin] [--nb-coh 15]
"""

import argparse
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal

sys.path.insert(0, '/home/shoshi/MITgcm_c69j/lab_sea12/mappings')
import plot_swot_spectra as pss

RUNSROOT = pss.RUNSROOT
HR = RUNSROOT + 'assim_swot_HR_10days_ctrlWyearly/'
MG = RUNSROOT + 'assim_swot_MG_10days_ctrlWyearly/'

GREEN, BLUE, CRIMSON = '#1baf7a', 'royalblue', 'crimson'

N_BINS_COH = 15      # log bins for the coherence panel; see docstring point 2
MIN_NEFF   = 3.0     # drop a coherence bin below this effective DOF

# Below this wavelength the HR model is dissipation-limited; shaded in both
# panels.  Measured from grid_hires_cleanbathy rather than assumed:  DYC is
# uniform at 9.266 km and DXC is 5.131 km (median) at the SWOT sample
# locations, and these tracks run at a median 69.5 deg from east -- nearly
# north-south -- so the spacing PROJECTED along the track is 8.864 km
# (p5 8.40, p95 9.07), dominated by DYC and almost double DXC.  50 km is
# 5.6 of those.  Note pss.MODEL_RES is 60 km and its comment reasons from the
# 4.65 km zonal spacing, which is the wrong metric for an along-track
# spectrum; left alone so the mg_20 figures do not shift.
MODEL_RES  = 50.0    # km

# (label, run dir, iteration, colour, linewidth, linestyle)
VARIANTS = {
    'requested': [
        ('Model iter 0  (both runs)',           MG,  0, CRIMSON, 1.5, ':'),
        ('HR-only  it0008  (m1qn3 iter 8)',     HR,  8, GREEN,   1.8, '--'),
        ('HR-only  it0020  (m1qn3 iter 18)',    HR, 20, GREEN,   2.2, '-'),
        ('Multi-grid  it0008  (rejected trial)', MG, 8, BLUE,    1.8, '--'),
    ],
    # lowest costfinal in each run, and in both cases also the last accepted
    # m1qn3 iterate:  HR it0020 = 2.0676e5 (iter 18),  MG it0020 = 3.3034e5
    # (iter 13).  Solid in both, since the iteration channel is not in use.
    'best': [
        ('Model iter 0  (both runs)',                       MG,  0, CRIMSON, 1.5, ':'),
        ('HR-only  it0020   (best, $f$ = 2.07$\\times$10$^5$)',   HR, 20, GREEN,   2.2, '-'),
        ('Multi-grid  it0020   (best, $f$ = 3.30$\\times$10$^5$)', MG, 20, BLUE,    2.2, '-'),
    ],
    'accepted': [
        ('Model iter 0  (both runs)',           MG,  0, CRIMSON, 1.5, ':'),
        ('HR-only  it0008  (m1qn3 iter 8)',     HR,  8, GREEN,   1.8, '--'),
        ('HR-only  it0020  (m1qn3 iter 18)',    HR, 20, GREEN,   2.2, '-'),
        ('Multi-grid  it0011  (m1qn3 iter 8)',  MG, 11, BLUE,    1.8, '--'),
        ('Multi-grid  it0020  (m1qn3 iter 13)', MG, 20, BLUE,    2.2, '-'),
    ],
    # mg_27.  Asymmetric on purpose -- the two iterations the comparison asked
    # for.  MG it0005 is m1qn3 iter 6 (f = 3.8716e5), not MG's best.
    'coherence': [
        ('Model iter 0',       MG,  0, CRIMSON, 1.6, ':'),
        ('Multi-grid  it0005', MG,  5, BLUE,    2.2, '-'),
        ('HR-only  it0020',    HR, 20, GREEN,   2.2, '-'),
    ],
}


# ========================== coherence machinery ==========================

def _seg_spec(x):
    """Tapered, detrended complex spectrum of one segment, mean bin dropped.

    No density normalization: every factor common to Sxx, Syy and Sxy cancels in
    the coherence ratio, so this deliberately does not match
    plot_swot_spectra.psd_on_master's 'density' scaling.
    """
    n = len(x)
    X = np.fft.rfft(signal.detrend(x) * signal.get_window('hann', n))[1:]
    return np.fft.rfftfreq(n, pss.DX)[1:], X


def coherence(segments, a, b, edges, nb):
    """Magnitude-squared coherence per bin, and the raw DOF count n_d.

    n_d is the number of (segment, native frequency) pairs that landed in the
    bin -- the count of independent complex products IF segments were
    independent.  See docstring point 3 for why they are not.
    """
    Sxx = np.zeros(nb)
    Syy = np.zeros(nb)
    Sxy = np.zeros(nb, dtype=complex)
    n_d = np.zeros(nb)
    for s in segments:
        f, X = _seg_spec(a[s])
        _, Y = _seg_spec(b[s])
        i = np.digitize(f, edges) - 1
        ok = (i >= 0) & (i < nb)
        i = i[ok]
        np.add.at(Sxx, i, (np.abs(X) ** 2)[ok])
        np.add.at(Syy, i, (np.abs(Y) ** 2)[ok])
        np.add.at(Sxy, i, (X * np.conj(Y))[ok])
        np.add.at(n_d, i, 1)
    coh2 = np.where(n_d > 0,
                    np.abs(Sxy) ** 2 / np.maximum(Sxx * Syy, 1e-300), np.nan)
    return coh2, n_d


def swath_redundancy(ds, segments):
    """n_swaths / n_segments -- the factor that turns n_d into n_eff.

    Segments from the same (cycle, pass_no) are the ten track_id lines of one
    swath, 10 km apart, so they carry nothing like ten independent
    realizations.  Collapsing each swath to one is pessimistic by construction;
    it brackets the truth from below together with the 1/n_d line.
    """
    key = np.stack([ds.cycle.values, ds.pass_no.values], axis=1)
    swaths = {tuple(key[s[0]]) for s in segments}
    return len(swaths) / float(len(segments)), len(swaths)


def main_coherence(outfile, curves, n_bins_psd, method_psd, nb_coh=N_BINS_COH):
    ds = xr.open_dataset(pss.OBS_FILE)
    nobs = ds.dims['iOBS']
    obs = ds.obs_val.values

    loaded = []
    for label, run, it, color, lw, ls in curves:
        mv, mm = pss.load_equi(run + 'run_adhi_it%04d/' % it, nobs)
        loaded.append((label, mv, mm, color, lw, ls))

    # every curve must sit on the SAME segments, so require all of them usable
    usable = np.isfinite(obs)
    for _, mv, mm, _, _, _ in loaded:
        usable &= (mm == 1) & np.isfinite(mv)
    segments = pss.build_segments(ds, usable)
    r_swath, n_swaths = swath_redundancy(ds, segments)
    print('usable samples: %d of %d' % (usable.sum(), nobs))
    print('gap-free segments with >= %d points: %d  (%d points, %d swaths)'
          % (pss.MIN_LEN, len(segments), sum(len(s) for s in segments), n_swaths))
    print('swath redundancy n_eff/n_d = %.3f' % r_swath)

    # ---- top panel: PSD on the mg_12b grid --------------------------------
    pss.set_grid(n_bins_psd, method_psd)
    f_psd = pss.master_frequencies
    psd_obs = pss.stack(segments, obs)
    psd = [(label, pss.stack(segments, mv), color, lw, ls)
           for label, mv, mm, color, lw, ls in loaded]

    # ---- bottom panel: coherence on its own coarse log-bin grid -----------
    pss.set_grid(nb_coh, 'logbin')
    f_coh, edges = pss.master_frequencies, pss._edges
    coh = []
    n_d_ref = None
    for label, mv, mm, color, lw, ls in loaded:
        c2, n_d = coherence(segments, obs, mv, edges, nb_coh)
        n_d_ref = n_d if n_d_ref is None else n_d_ref
        coh.append((label, c2, color, lw, ls))
    n_eff = n_d_ref * r_swath
    keep = n_eff >= MIN_NEFF
    floor_opt = np.where(n_d_ref > 0, 1.0 / np.maximum(n_d_ref, 1), np.nan)
    floor_eff = np.where(keep, 1.0 / np.maximum(n_eff, 1e-30), np.nan)
    coh = [(lb, np.where(keep, c2, np.nan), c, lw, ls) for lb, c2, c, lw, ls in coh]

    lams = (1000, 500, 200, 150, 100, 70, 50, 35, 25)
    print('\n=== coherence^2 with SWOT, %d log bins ===' % nb_coh)
    hdr = '%-20s' % 'lam(km)' + ''.join('%-9d' % l for l in lams)
    print(hdr)

    def row(name, arr, fmt='%-9.3f'):
        i = [np.argmin(np.abs(1 / f_coh - l)) for l in lams]
        print('%-20s' % name + ''.join(fmt % arr[j] for j in i))

    for label, c2, _, _, _ in coh:
        row(label.split('(')[0].strip(), c2)
    row('n_d', n_d_ref, '%-9.0f')
    row('n_eff', n_eff, '%-9.1f')
    row('floor 1/n_eff', floor_eff)

    # The number worth quoting: the shortest scale down to which a run still
    # holds coherence distinguishable from zero.  Take the block that is
    # CONTIGUOUS FROM THE LONG-WAVELENGTH END, so an isolated bin poking above
    # the floor further down cannot extend the claim.
    print('\ncoherence above the 1/n_eff floor, contiguously, down to:')
    for label, c2, _, _, _ in coh:
        above = np.where(np.isfinite(c2) & (c2 > floor_eff))[0]
        gap = np.where(np.diff(above) != 1)[0]
        blk = above[:gap[0] + 1] if len(gap) else above
        print('  %-40s %s' % (label.split('(')[0].strip(),
                              ('%.0f km' % (1 / f_coh[blk[-1]])) if len(blk)
                              else 'never'))

    # ------------------------------------------------------------------ plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10.2), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1.55, 1],
                                                    hspace=0.08))
    ms_psd = 3 if n_bins_psd <= 60 else 0

    mean_ds, ci_ds = psd_obs[0], psd_obs[1]
    ax1.loglog(f_psd, mean_ds, color='black', lw=2.2, marker='o', ms=ms_psd,
               label='SWOT observations')
    ax1.fill_between(f_psd, np.maximum(mean_ds - ci_ds, 1e-10), mean_ds + ci_ds,
                     color='black', alpha=0.15)
    for label, (mean, ci, _, _), color, lw, ls in psd:
        ax1.loglog(f_psd, mean, color=color, lw=lw, ls=ls, marker='o',
                   ms=ms_psd, label=label)
        ax1.fill_between(f_psd, np.maximum(mean - ci, 1e-10), mean + ci,
                         color=color, alpha=0.15)
    ax1.axvspan(1 / MODEL_RES, pss.F_NYQ, color='0.85', alpha=0.5, zorder=0)
    ax1.text(0.033, 0.055, 'below model\neffective res.', fontsize=8,
             color='0.35', ha='center', transform=ax1.get_xaxis_transform())
    ax1.set_ylabel('PSD (m$^2$ / cycles/km)')
    ax1.set_ylim(1e-6, 2e2)
    ax1.grid(True, which='both', ls='--', alpha=0.4)
    ax1.legend(loc='lower left', fontsize=9)
    ax1.set_title('SWOT vs model along-track SSH, 10-day window, '
                  'ctrlWyearly runs', pad=24, fontsize=11.5)

    axt = ax1.secondary_xaxis('top', functions=(lambda x: 1 / x, lambda x: 1 / x))
    axt.set_xlabel('Wavelength (km)')
    axt.set_xticks([2000, 1000, 500, 200, 100, 50, 30, 20])
    axt.set_xticklabels(['2000', '1000', '500', '200', '100', '50', '30', '20'])

    for label, c2, color, lw, ls in coh:
        ax2.semilogx(f_coh, c2, color=color, lw=lw, ls=ls, marker='o', ms=4,
                     zorder=3, label=label)
    ax2.axvspan(1 / MODEL_RES, pss.F_NYQ, color='0.85', alpha=0.5, zorder=0)
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel('coherence$^2$ with SWOT')
    ax2.set_xlabel('Wavenumber (cycles / km)')
    ax2.grid(True, which='both', ls='--', alpha=0.4)
    ax2.legend(loc='upper right', fontsize=8.5, framealpha=0.92)
    ax2.text(0.033, 0.06, 'below model\neffective res.', fontsize=8,
             color='0.35', ha='center', transform=ax2.get_xaxis_transform())

    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    print('\nwrote %s' % outfile)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('outfile', nargs='?', default=None)
    ap.add_argument('--variant', choices=sorted(VARIANTS), default='requested')
    ap.add_argument('-n', '--n-bins', type=int, default=120)
    ap.add_argument('-m', '--method', choices=['logbin', 'interp'], default='interp')
    ap.add_argument('--nb-coh', type=int, default=N_BINS_COH,
                    help='log bins for the coherence panel (--variant coherence)')
    ap.add_argument('--no-ratio', action='store_true')
    a = ap.parse_args()
    out = a.outfile or {
        'requested': 'mg_20b_swot_spectra_HR_vs_MG_ctrlWyearly.png',
        'accepted':  'mg_20c_swot_spectra_HR_vs_MG_ctrlWyearly_accepted.png',
        'best':      'mg_20d_swot_spectra_HR_vs_MG_ctrlWyearly_best.png',
        'coherence': 'mg_27_swot_coherence_HR_vs_MG_ctrlWyearly.png',
    }[a.variant]
    if a.variant == 'coherence':
        main_coherence(out, VARIANTS['coherence'], a.n_bins, a.method, a.nb_coh)
    else:
        pss.main(out, a.n_bins, a.method, VARIANTS[a.variant], not a.no_ratio)
