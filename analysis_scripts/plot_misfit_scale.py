#!/usr/bin/env python3
"""
At what spatial scale does the SWOT assimilation actually help, and what sets
that scale?

Three configurations that isolate the question:

    LR   1/3 deg forward  + 1/3 deg adjoint
    MG   1/12 deg forward + 1/3 deg adjoint      <- multi-grid
    HR   1/12 deg forward + 1/12 deg adjoint

MG and HR share a forward model and are bit-identical at iteration 0, so they
differ in exactly one thing: adjoint resolution.  That is the controlled
experiment for "does the coarse adjoint cap the scale of the correction?"
LR additionally changes the forward model, so it answers a different question.

Six panels:
  (a) misfit PSD per configuration.
  (b) misfit PSD ratio to iteration 0.  Below 1 the optimizer removed variance,
      above 1 it added some; the crossing is where assimilation stops helping.
  (c) coherence^2 between SWOT and the model -- PSD says whether the model has
      the right AMOUNT of variance at a scale, coherence says whether it is in
      the right PLACE.  Variance added out of phase raises the misfit.
  (d) increment spectrum, PSD(mod_it5 - mod_it0), normalised at 400 km.  THE
      test of the hypothesis: if adjoint resolution caps the correction scale,
      MG must fall off faster than HR here.
  (e) model SSH PSD / observed -- how much variance the model even has.
  (f) misfit rms and variance-weighted centroid wavelength per iteration,
      lam_c = 1 / (int f*S df / int S df).

Coherence is biased high by ~1/n_eff, and n_eff is NOT the segment count:
track_id 0-9 are ten parallel lines 10 km apart in one swath, over 3 repeat
cycles.  Treat anything below ~0.05 as indistinguishable from zero.

Usage:  python3 plot_misfit_scale.py [outfile.png] [--weights v1|v2]
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
from plot_swot_spectra import (load_equi, build_segments, stack, set_grid,
                               OBS_FILE, RUNSROOT, DX)
import plot_swot_spectra as P

# ============================== configuration ==============================
def runs_for(weights):
    """weights='v2' uses the August-2026 regenerated control priors."""
    sfx = '_v2' if weights == 'v2' else ''
    return {
        'LR': (RUNSROOT + 'assim_swot_LR_10days/', 'run_adlo_it%04d/'),
        'MG': (RUNSROOT + 'assim_swot_MG_10days_xxmult16_newCtrlW%s/' % sfx,
               'run_adhi_it%04d/'),
        'HR': (RUNSROOT + 'assim_swot_HR_10days_newCtrlW%s/' % sfx,
               'run_adhi_it%04d/'),
    }

# Per-run reference iteration for panels (a), (c), (d), (e).  m1qn3 evaluates
# the cost at trial points too, so it#### is not automatically an iterate --
# these are the last ACCEPTED ones (MG v1 it5, the old shared IT_REF, is a
# rejected trial).  Verify with plot_mg_diagnostics.classify_iterates.
def itref_for(weights):
    return ({'LR': 5, 'MG': 5, 'HR': 10} if weights == 'v2'
            else {'LR': 5, 'MG': 8, 'HR': 17})

MAXIT    = {'LR': 10, 'MG': 10, 'HR': 10}          # panel (f) trajectory, per run
SHOW_ITS = [1, 2, 3, 5, 10]                        # panel (b)
N_BINS   = 40
LAM_MAX  = 700.                                    # panels b-e: >90 segments contribute

COL   = {'LR': '#eda100', 'MG': '#4169e1', 'HR': '#1baf7a'}
LS    = {'LR': '-.',      'MG': '-',       'HR': '-'}
CIT0  = '#dc143c'
RAMP  = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#0d366b']   # ordinal blue
LR_DX, LR_DY = 20.7, 37.1                          # km at 56N, 1/3 deg
NOISE = 2 * (0.02 * 0.208) ** 2 * DX               # KaRIn white-noise PSD
# ===========================================================================

set_grid(N_BINS, 'logbin')
F, EDGES, NB = P.master_frequencies, P._edges, N_BINS


def _spec(x):
    N = len(x)
    X = np.fft.rfft(signal.detrend(x) * signal.get_window('hann', N))[1:]
    return np.fft.rfftfreq(N, DX)[1:], X


def coherence(segments, a, b):
    """Magnitude-squared coherence.  Cross- and auto-spectra must be summed
    BEFORE the ratio is formed -- per-segment coherence is identically 1."""
    Sxx = np.zeros(NB); Syy = np.zeros(NB)
    Sxy = np.zeros(NB, dtype=complex); cnt = np.zeros(NB)
    for s in segments:
        f, X = _spec(a[s]); _, Y = _spec(b[s])
        i = np.digitize(f, EDGES) - 1
        ok = (i >= 0) & (i < NB); i = i[ok]
        np.add.at(Sxx, i, (np.abs(X) ** 2)[ok])
        np.add.at(Syy, i, (np.abs(Y) ** 2)[ok])
        np.add.at(Sxy, i, (X * np.conj(Y))[ok])
        np.add.at(cnt, i, 1)
    return np.where(cnt > 0, np.abs(Sxy) ** 2 / np.maximum(Sxx * Syy, 1e-300), np.nan)


def variance(S):
    m = np.isfinite(S)
    return np.trapz(S[m], F[m])


def centroid_lambda(S):
    m = np.isfinite(S)
    return 1.0 / (np.trapz(F[m] * S[m], F[m]) / np.trapz(S[m], F[m]))


def frac_below(S, lam):
    m = np.isfinite(S)
    hi = m & (F > 1.0 / lam)
    return np.trapz(S[hi], F[hi]) / np.trapz(S[m], F[m])


def crossing(ratio, target=1.0):
    """Longest wavelength beyond which the ratio stays above `target`.

    Scanned from the Nyquist end down, so an isolated wiggle at low wavenumber
    -- where only a few dozen segments contribute -- cannot masquerade as the
    crossover.  NaN if the ratio never exceeds `target` at the small-scale end.
    """
    m = np.isfinite(ratio)
    f, r = F[m], ratio[m]
    if r[-1] <= target:
        return np.nan
    for j in range(len(f) - 1, 0, -1):
        if r[j - 1] <= target:
            w = (target - r[j - 1]) / (r[j] - r[j - 1])
            return 1.0 / np.exp(np.log(f[j - 1]) + w * (np.log(f[j]) - np.log(f[j - 1])))
    return np.nan


def at(S, lam):
    return S[np.argmin(np.abs(1 / F - lam))]


def main(outfile='misfit_scale.png', weights='v2'):
    RUNS = runs_for(weights)
    IT_REF = itref_for(weights)
    ITS = {t: sorted(set(range(0, MAXIT[t] + 1)) | {IT_REF[t]}) for t in RUNS}
    ds = xr.open_dataset(OBS_FILE)
    nobs = ds.dims['iOBS']
    obs = ds.obs_val.values

    mod = {}
    usable = np.isfinite(obs)
    for tag, (root, pat) in RUNS.items():
        for it in ITS[tag]:
            mv, mm = load_equi(root + pat % it, nobs)
            mod[(tag, it)] = mv
            usable &= (mm == 1) & np.isfinite(mv)
    segments = build_segments(ds, usable)
    print('control weights: %s' % weights)
    print('%d usable samples, %d segments, %d points'
          % (usable.sum(), len(segments), sum(len(s) for s in segments)))

    S_obs = stack(segments, obs)[0]
    mis = {k: stack(segments, obs - v)[0] for k, v in mod.items()}
    inc = {t: stack(segments, mod[(t, IT_REF[t])] - mod[(t, 0)])[0] for t in RUNS}
    smod = {(t, i): stack(segments, mod[(t, i)])[0]
            for t in RUNS for i in (0, IT_REF[t])}
    g2 = {(t, i): coherence(segments, obs, mod[(t, i)])
          for t in RUNS for i in (0, IT_REF[t])}

    lams = [400, 300, 200, 150, 100, 70, 50, 35, 25]
    hdr = '%-10s ' % 'lam(km)' + ''.join('%-8d' % l for l in lams)

    print('\n=== misfit PSD ratio to iteration 0 ===')
    print(hdr)
    for t in ('LR', 'MG', 'HR'):
        for it in sorted({IT_REF[t], MAXIT[t]}):
            r = mis[(t, it)] / mis[(t, 0)]
            print('%-10s ' % ('%s it%d' % (t, it)) + ''.join('%-8.2f' % at(r, l) for l in lams))

    print('\n=== increment PSD(it - it0), normalised at 400 km ===')
    print(hdr)
    for t in ('LR', 'MG', 'HR'):
        print('%-10s ' % (t + ' inc') + ''.join('%-8.3f' % (at(inc[t], l) / at(inc[t], 400)) for l in lams))

    print('\n=== coherence^2 with SWOT  (< ~0.05 = indistinguishable from zero) ===')
    print(hdr)
    for t in ('LR', 'MG', 'HR'):
        for it in (0, IT_REF[t]):
            print('%-10s ' % ('%s it%d' % (t, it)) + ''.join('%-8.2f' % at(g2[(t, it)], l) for l in lams))

    print('\n=== model SSH PSD / observed ===')
    print(hdr)
    for t in ('LR', 'MG', 'HR'):
        for it in (0, IT_REF[t]):
            print('%-10s ' % ('%s it%d' % (t, it)) +
                  ''.join('%-8.2f' % (at(smod[(t, it)], l) / at(S_obs, l)) for l in lams))

    print('\n=== summary ===')
    print('%-4s %-6s %-9s %-9s %-16s %-11s %s'
          % ('run', 'itref', 'rms it0', 'rms itref', 'lam_c it0->itref',
             'frac<100km', 'ratio=1 at (km)'))
    for t in ('LR', 'MG', 'HR'):
        r = IT_REF[t]
        x = crossing(mis[(t, r)] / mis[(t, 0)])
        print('%-4s %-6d %-9.4f %-9.4f %-16s %-11s %s'
              % (t, r, np.sqrt(variance(mis[(t, 0)])), np.sqrt(variance(mis[(t, r)])),
                 '%.0f -> %.0f' % (centroid_lambda(mis[(t, 0)]),
                                   centroid_lambda(mis[(t, r)])),
                 '%.3f -> %.3f' % (frac_below(mis[(t, 0)], 100.),
                                   frac_below(mis[(t, r)], 100.)),
                 '%.0f' % x if np.isfinite(x) else 'never above 1'))

    # ------------------------------------------------------------------ plot
    fig, ax = plt.subplots(2, 3, figsize=(19, 10))
    ticks = [1000, 500, 200, 100, 50, 25]

    def deco(a, lam_max=1600., top=True):
        a.set_xscale('log'); a.set_xlim(1.0 / lam_max, 5.5e-2)
        a.grid(True, which='both', ls='--', alpha=0.35)
        a.set_xlabel('Wavenumber (cycles / km)')
        if top:
            t = a.secondary_xaxis('top', functions=(lambda x: 1/x, lambda x: 1/x))
            t.set_xticks(ticks); t.set_xticklabels([str(v) for v in ticks])
            t.set_xlabel('Wavelength (km)')

    # (a) misfit spectrum
    a = ax[0, 0]
    a.loglog(F, S_obs, color='0.35', lw=1.3, ls=':', label='SWOT (reference)')
    a.loglog(F, mis[('MG', 0)], color=CIT0, lw=2.0, ls='--', label='misfit, iter 0 (MG/HR)')
    for t in ('LR', 'MG', 'HR'):
        a.loglog(F, mis[(t, IT_REF[t])], color=COL[t], lw=2.0, ls=LS[t],
                 label='misfit, %s iter %d' % (t, IT_REF[t]))
    a.axhline(NOISE, color='0.5', ls=':', lw=1.3)
    a.set_ylabel('misfit PSD (m$^2$ / cycles/km)')
    a.set_title('(a) misfit spectrum', loc='left', fontsize=11)
    a.legend(fontsize=8, loc='lower left'); deco(a)

    # (b) ratio to iteration 0
    a = ax[0, 1]
    for c, it in zip(RAMP, SHOW_ITS):
        a.semilogx(F, mis[('MG', it)] / mis[('MG', 0)], color=c, lw=1.8,
                   label='MG iter %d' % it)
    for t in ('HR', 'LR'):
        a.semilogx(F, mis[(t, IT_REF[t])] / mis[(t, 0)], color=COL[t], lw=1.8,
                   ls=LS[t] if t == 'LR' else '--', label='%s iter %d' % (t, IT_REF[t]))
    a.axhline(1.0, color='0.25', lw=1.1)
    x = crossing(mis[('MG', IT_REF['MG'])] / mis[('MG', 0)])
    a.axvline(1 / x, color='0.45', ls=':', lw=1.4)
    a.annotate('MG stops\nhelping: %.0f km' % x, xy=(1/x, 0.52), xytext=(1/x/3.0, 0.40),
               fontsize=8.5, color='0.3', ha='center',
               arrowprops=dict(arrowstyle='->', color='0.5', lw=0.9))
    a.axvspan(1/(2*LR_DY), 1/(2*LR_DX), color='0.85', alpha=0.6, zorder=0)
    a.text(1/(2*LR_DX)*0.97, 0.36, '2$\\Delta$ LR grid', fontsize=7.5, color='0.35',
           ha='right', rotation=90, va='bottom')
    a.set_ylim(0.3, 1.85); a.set_ylabel('misfit PSD / iteration 0')
    a.set_title('(b) where the optimizer helps (<1) and hurts (>1)', loc='left', fontsize=11)
    a.legend(fontsize=8, loc='upper left', ncol=2); deco(a, LAM_MAX)

    # (c) coherence
    a = ax[0, 2]
    a.semilogx(F, g2[('MG', 0)], color=CIT0, lw=2.0, ls='--', marker='o', ms=3,
               label='iter 0 (MG/HR)')
    for t in ('LR', 'MG', 'HR'):
        a.semilogx(F, g2[(t, IT_REF[t])], color=COL[t], lw=2.0, ls=LS[t], marker='o', ms=3,
                   label='%s iter %d' % (t, IT_REF[t]))
    a.axhspan(0, 0.05, color='0.85', alpha=0.7, zorder=0)
    a.text(1.6e-3, 0.015, 'indistinguishable from zero', fontsize=8, color='0.35')
    a.set_ylim(0, 0.85); a.set_ylabel('coherence$^2$ (SWOT vs model)')
    a.set_title('(c) is the variance in the right PLACE?', loc='left', fontsize=11)
    a.legend(fontsize=8); deco(a, LAM_MAX)

    # (d) increment spectrum -- the hypothesis test
    a = ax[1, 0]
    for t in ('LR', 'MG', 'HR'):
        a.loglog(F, inc[t] / at(inc[t], 400), color=COL[t], lw=2.0, ls=LS[t],
                 marker='o', ms=3, label='%s' % t)
    a.set_ylim(1e-4, 3)
    a.set_ylabel('increment PSD, normalised at 400 km')
    a.set_title('(d) MG vs HR difference has no consistent sign $\\Rightarrow$ '
                'not an adjoint-resolution limit', loc='left', fontsize=10)
    a.legend(fontsize=8.5); deco(a, LAM_MAX)

    # (e) model variance
    a = ax[1, 1]
    a.axhline(1.0, color='0.25', lw=1.1)
    for t in ('LR', 'MG', 'HR'):
        a.loglog(F, smod[(t, IT_REF[t])] / S_obs, color=COL[t], lw=2.0, ls=LS[t],
                 label='%s iter %d' % (t, IT_REF[t]))
    a.loglog(F, smod[('MG', 0)] / S_obs, color=CIT0, lw=2.0, ls='--',
             label='iter 0 (MG/HR)')
    a.loglog(F, smod[('LR', 0)] / S_obs, color=COL['LR'], lw=1.2, ls=':',
             label='iter 0 (LR)')
    a.set_ylim(2e-3, 3); a.set_ylabel('model SSH PSD / observed')
    a.set_title('(e) how much variance the model even has', loc='left', fontsize=11)
    a.legend(fontsize=8, loc='lower left'); deco(a, LAM_MAX)

    # (f) scalars per iteration
    a = ax[1, 2]
    for t in ('LR', 'MG', 'HR'):
        tr = [i for i in ITS[t] if i <= MAXIT[t]]
        a.plot(tr, [centroid_lambda(mis[(t, i)]) for i in tr], color=COL[t],
               lw=2.0, ls=LS[t], marker='o', ms=4, label='%s  $\\lambda_c$' % t)
    a.set_xlabel('outer iteration'); a.set_ylabel('misfit centroid wavelength (km)')
    a.grid(True, ls='--', alpha=0.35)
    a2 = a.twinx()
    for t in ('LR', 'MG', 'HR'):
        tr = [i for i in ITS[t] if i <= MAXIT[t]]
        a2.plot(tr, [np.sqrt(variance(mis[(t, i)])) for i in tr], color=COL[t],
                lw=1.3, ls=':', marker='s', ms=3)
    a2.set_ylabel('misfit rms (m)   [dotted]')
    a.set_title('(f) misfit gets smaller AND smaller-scale', loc='left', fontsize=11)
    a.legend(fontsize=8.5, loc='upper right')

    fig.suptitle('What sets the spatial scale of SWOT assimilation skill?   '
                 '10-day window, %d along-track segments, control weights %s'
                 % (len(segments), weights), fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(outfile, dpi=150)
    print('\nwrote %s' % outfile)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('outfile', nargs='?', default='misfit_scale.png')
    ap.add_argument('-w', '--weights', choices=['v1', 'v2'], default='v2')
    a = ap.parse_args()
    main(a.outfile, a.weights)
