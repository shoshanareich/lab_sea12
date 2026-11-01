#!/usr/bin/env python3
"""
Along-track wavenumber spectra of SWOT SSH vs the model, over outer iterations.

Compares the observations against the model equivalents that pkg/obsfit writes
at the SWOT sample locations, for the multi-grid run and the single-resolution
HR run.  Iteration 0 is bit-identical between the two (verified: max |diff| = 0
over all 27960 usable samples), so it is drawn once as the shared baseline.

Method
------
1. Group samples by (cycle, pass_no, track_id).  track_id 0-9 is the CROSS-track
   pixel and line_id 0-270 is the along-track index, so this key selects one
   along-track line -- the right unit for a 1-D spectrum.
2. Split each line into GAP-FREE segments on a discontinuity in line_id.  A
   dropped sample (land, quality flag, outside the run window) leaves two points
   more than 10 km apart, and periodogram assumes uniform spacing.  58% of the
   lines here contain at least one gap, up to 650 km.  Ignoring them inflates the
   observed PSD by 2.7x at 50 km and 10.5x at 25 km.
3. Keep segments of >= MIN_LEN points.  Short segments carry no low-wavenumber
   information and the Hann taper eats their ends, so they only add chi2_2 noise
   to the high-wavenumber bins.
4. Periodogram each segment (detrend -> Hann -> density), then AVERAGE the native
   values within each log-spaced bin.  Averaging m native values gives chi2_2m
   rather than chi2_2, and the bins are independent -- interpolating onto the
   master grid instead leaves adjacent bins correlated and, at low wavenumber,
   invents values where the segment has no Fourier content.
5. Mean and 95% interval across segments, per bin.

Reading the result
------------------
* obs_val is SSH, not SSHA -- the preprocessing does ssh = ssha_filtered + mdt.
  So lambda > ~200 km is dominated by mean dynamic topography, and there you are
  comparing the model mean state against the MDT, not eddy variance.
* Nyquist here is 20 km, not SWOT's native 2 km: the obs file was Gaussian
  smoothed (sigma = 2 px ~ 4 km, truncated at +/-2 px) then decimated [5:-5:5].
  That filter leaves ~20% of the power at the new Nyquist, so the top couple of
  bins carry aliased energy.
* The model cannot hold energy below ~6 dx_HR (~30-60 km, latitude dependent).
  The roll-off there is numerical dissipation, not a misfit.
* N in the CI is the segment count, but track_id 0-9 are ten parallel lines only
  10 km apart in one swath, over 3 repeat cycles -- the true DOF is far smaller,
  so the shading is a lower bound on the uncertainty.

Usage:  python3 plot_swot_spectra.py [outfile.png]
"""

import argparse
import glob
import sys

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal

# ============================== configuration ==============================
RUNSROOT = '/scratch/shoshi/labsea_MG_12/'
OBS_FILE = ('/home/shoshi/MITgcm_c69j/lab_sea12/input_binaries_hires/'
            'swot_obsfit_cycles_9thru11_labsea_L3v3.nc')
OBS_NAME = 'swot_obsfit_cycles_9thru11_labsea_L3v3'

# (label, run dir, iteration, colour, linewidth, linestyle)
def models_for(weights='v2', it_mg=5, it_hr=5):
    """weights='v2' selects the August-2026 regenerated control priors.

    Iteration numbers are per run because m1qn3 evaluates the cost at trial
    points as well as at iterates, so it#### does not mean the same thing in
    two runs -- check plot_mg_diagnostics.classify_iterates before comparing.
    """
    sfx = '_v2' if weights == 'v2' else ''
    mg = RUNSROOT + 'assim_swot_MG_10days_xxmult16_newCtrlW%s/' % sfx
    hr = RUNSROOT + 'assim_swot_HR_10days_newCtrlW%s/' % sfx
    return [
        ('Model iter 0 (both runs)',        mg, 0,     'crimson',   1.4, '--'),
        ('Multi-grid iter %d' % it_mg,      mg, it_mg, 'royalblue', 1.8, '-'),
        ('HR-only iter %d' % it_hr,         hr, it_hr, '#1baf7a',   1.8, '-'),
    ]

MODELS = models_for()

# hue carries the configuration, line style the control-prior version, so both
# comparisons (MG vs HR, v1 vs v2) read off separate channels rather than four
# competing hues.  crimson/royalblue/#1baf7a validate all-pairs in light mode.
_CFG = {'MG': ('assim_swot_MG_10days_xxmult16_newCtrlW', 'royalblue', 'Multi-grid'),
        'HR': ('assim_swot_HR_10days_newCtrlW',          '#1baf7a',   'HR-only')}
_VER = {'v1': ('--', 1.7), 'v2': ('-', 2.0)}


def models_from_spec(spec):
    """Parse 'MG:v1:8,MG:v2:5,HR:v1:17,HR:v2:10' into a MODELS list.

    Iteration 0 is bit-identical across all four runs (verified: max |diff| = 0
    over all 27960 usable samples), so it is drawn once from the first entry.
    """
    out, first = [], None
    for item in spec.split(','):
        cfg, ver, it = item.split(':')
        base, color, name = _CFG[cfg]
        ls, lw = _VER[ver]
        run = RUNSROOT + base + ('_v2' if ver == 'v2' else '') + '/'
        first = first or run
        out.append(('%s %s  iter %s' % (name, ver, it), run, int(it), color, lw, ls))
    return [('Model iter 0 (all runs)', first, 0, 'crimson', 1.4, ':')] + out

DX        = 10.0                  # along-track sample spacing, km (measured 9.9855)
MIN_LEN   = 30                    # points; shorter segments are dropped
F_MIN     = 0.0005                # cycles/km -> 2000 km
N_BINS    = 40                    # log bins; 120 leaves most bins empty below 200 km
SIGMA_OBS = 0.02 * 0.208          # KaRIn noise after the sigma=2px Gaussian, m
MODEL_RES = 60.0                  # km; below this the model is dissipation-limited
# ===========================================================================

F_NYQ  = 1.0 / (2.0 * DX)
METHOD = 'logbin'                 # 'logbin' | 'interp' (the original behaviour)

master_frequencies = None
_edges = None


def set_grid(n_bins=N_BINS, method=METHOD):
    """(Re)build the master frequency grid and its bin edges."""
    global master_frequencies, _edges, N_BINS, METHOD
    N_BINS, METHOD = n_bins, method
    master_frequencies = np.logspace(np.log10(F_MIN), np.log10(F_NYQ), N_BINS)
    _mid = np.sqrt(master_frequencies[:-1] * master_frequencies[1:])
    _edges = np.r_[master_frequencies[0] ** 2 / _mid[0], _mid,
                   master_frequencies[-1] ** 2 / _mid[-1]]


set_grid()


def load_equi(run_dir, nobs, obs_name=OBS_NAME):
    """Gather the per-tile .equi.nc files into arrays aligned with the obs file.

    obsfit writes one file per MPI tile holding only that tile's samples, keyed
    by the 1-based global sample index.  mod_mask stays 0 where a sample was
    never written -- outside the run window, or on a dry/edge cell.  Note the
    tile files declare mod_val:_FillValue = 0., so initialise to NaN rather than
    letting a masked sample read back as a plausible-looking zero.
    """
    mod_val = np.full(nobs, np.nan)
    mod_mask = np.zeros(nobs)
    files = sorted(glob.glob(run_dir + obs_name + '.*.equi.nc'))
    if not files:
        raise IOError('no %s.*.equi.nc under %s' % (obs_name, run_dir))
    for fn in files:
        d = xr.open_dataset(fn)
        i = d.sample_ind_glob.values.astype(int) - 1
        ok = (i >= 0) & (i < nobs)
        mod_val[i[ok]] = d.mod_val.values[ok]
        mod_mask[i[ok]] = d.mod_mask.values[ok]
        d.close()
    return mod_val, mod_mask


def build_segments(ds, usable, min_len=MIN_LEN):
    """Gap-free along-track segments of usable samples, >= min_len points."""
    groups = pd.DataFrame({'cycle': ds.cycle.values,
                           'pass_no': ds.pass_no.values,
                           'track_id': ds.track_id.values}
                          ).groupby(['cycle', 'pass_no', 'track_id']).groups
    line_id = ds.line_id.values
    segments = []
    for combo, indices in groups.items():
        idx = np.asarray(indices)
        idx = idx[usable[idx]]
        if len(idx) < 2:
            continue
        for seg in np.split(idx, np.where(np.diff(line_id[idx]) != 1)[0] + 1):
            if len(seg) >= min_len:
                segments.append(seg)
    return segments


def psd_on_master(x):
    """Periodogram of one segment, averaged within each log bin.

    A bin with no native frequency in it stays NaN, which at the low-wavenumber
    end replaces an explicit  f < 1/(N*dx)  mask: a segment simply has no
    Fourier content below its own fundamental.
    """
    freqs, native = signal.periodogram(signal.detrend(x), fs=1.0 / DX,
                                       window='hann', scaling='density')
    freqs, native = freqs[1:], native[1:]              # drop the mean

    if METHOD == 'interp':
        # Original behaviour: geometric interpolation onto the master grid, with
        # an explicit mask below the segment's fundamental.  Adjacent bins are
        # then interpolated from the same native pair, so they are correlated --
        # the curve looks smoother, and the CI narrower, than the data supports.
        y = np.exp(np.interp(master_frequencies, freqs, np.log(native)))
        y[master_frequencies < 1.0 / (len(x) * DX)] = np.nan
        return y, (~np.isnan(y)).astype(int)

    b = np.digitize(freqs, _edges) - 1
    ok = (b >= 0) & (b < N_BINS)
    b, native = b[ok], native[ok]
    count = np.bincount(b, minlength=N_BINS)
    total = np.bincount(b, weights=native, minlength=N_BINS)
    return np.where(count > 0, total / np.maximum(count, 1), np.nan), count


def stack(segments, values):
    """Mean PSD over segments, plus the 95% interval and the occupancy counts."""
    out = [psd_on_master(values[s]) for s in segments]
    M = np.array([o[0] for o in out])              # (n_segments, n_bins)
    ndof = np.array([o[1] for o in out]).sum(0)    # native values averaged, per bin
    mean = np.nanmean(M, axis=0)
    std = np.nanstd(M, axis=0)
    nseg = np.sum(~np.isnan(M), axis=0)
    return mean, 1.96 * std / np.sqrt(np.maximum(nseg, 1)), nseg, ndof


def main(outfile='swot_spectra.png', n_bins=N_BINS, method=METHOD, models=None,
         ratio=False):
    set_grid(n_bins, method)
    models = MODELS if models is None else models
    ds = xr.open_dataset(OBS_FILE)
    nobs = ds.dims['iOBS']
    obs = ds.obs_val.values

    loaded = []
    for label, run, it, color, lw, ls in models:
        mv, mm = load_equi(run + 'run_adhi_it%04d/' % it, nobs)
        loaded.append((label, mv, mm, color, lw, ls))

    # every curve must sit on the SAME segments, so require all of them usable
    usable = np.isfinite(obs)
    for _, mv, mm, _, _, _ in loaded:
        usable &= (mm == 1) & np.isfinite(mv)
    print('usable samples: %d of %d' % (usable.sum(), nobs))

    segments = build_segments(ds, usable)
    npts = sum(len(s) for s in segments)
    print('gap-free segments with >= %d points: %d  (%d points)'
          % (MIN_LEN, len(segments), npts))
    print('master grid: %d bins, %.1f%% per step, method=%s'
          % (N_BINS, 100 * (F_NYQ / F_MIN) ** (1.0 / (N_BINS - 1)) - 100, METHOD))

    mean_ds, ci_ds, nseg, ndof = stack(segments, obs)
    curves = [(label, stack(segments, mv), color, lw, ls)
              for label, mv, mm, color, lw, ls in loaded]

    hdr = '%-9s %-6s %-7s %-11s' % ('lam(km)', 'nseg', 'native', 'obs PSD')
    hdr += ''.join('%-14s' % (c[0].split('(')[0].strip() + '/obs') for c in curves)
    print('\n' + hdr)
    for lam in (1000, 500, 200, 150, 100, 70, 50, 35, 25):
        i = np.argmin(np.abs(1 / master_frequencies - lam))
        row = '%-9.0f %-6d %-7d %-11.3e' % (1 / master_frequencies[i],
                                            nseg[i], ndof[i], mean_ds[i])
        row += ''.join('%-14.2f' % (c[1][0][i] / mean_ds[i]) for c in curves)
        print(row)

    # ------------------------------------------------------------------ plot
    if ratio:
        fig, (ax1, axr) = plt.subplots(2, 1, figsize=(10, 9.6), sharex=True,
                                       gridspec_kw=dict(height_ratios=[2.1, 1],
                                                        hspace=0.07))
    else:
        fig, ax1 = plt.subplots(figsize=(9.5, 6.8))
        axr = None

    ms = 3 if N_BINS <= 60 else 0
    ax1.loglog(master_frequencies, mean_ds, color='black', lw=2.2,
               marker='o', ms=ms, label='SWOT observations')
    ax1.fill_between(master_frequencies, np.maximum(mean_ds - ci_ds, 1e-10),
                     mean_ds + ci_ds, color='black', alpha=0.15)
    for label, (mean, ci, _, _), color, lw, ls in curves:
        ax1.loglog(master_frequencies, mean, color=color, lw=lw, ls=ls,
                   marker='o', ms=ms, label=label)
        ax1.fill_between(master_frequencies, np.maximum(mean - ci, 1e-10),
                         mean + ci, color=color, alpha=0.15)

    ax1.axhline(2 * SIGMA_OBS ** 2 * DX, color='0.45', ls=':', lw=1.5,
                label='SWOT noise floor (%.0f mm)' % (1000 * SIGMA_OBS))
    ax1.axvspan(1 / MODEL_RES, F_NYQ, color='0.85', alpha=0.5, zorder=0)
    ax1.text(0.033, 0.055, 'below model\neffective res.', fontsize=8, color='0.35',
             ha='center', transform=ax1.get_xaxis_transform())

    ax1.set_xlabel('Wavenumber (cycles / km)')
    ax1.set_ylabel('PSD (m$^2$ / cycles/km)')
    ax1.set_title('SWOT vs model along-track SSH spectra, 10-day window\n'
                  '%d gap-free segments, %d bins, %s'
                  % (len(segments), N_BINS,
                     'log-bin averaged' if METHOD == 'logbin' else 'np.interp onto master grid'),
                  pad=22)
    ax1.grid(True, which='both', ls='--', alpha=0.4)
    ax1.legend(loc='lower left', fontsize=9)
    ax1.set_ylim(1e-6, 2e2)

    ax2 = ax1.secondary_xaxis('top', functions=(lambda x: 1 / x, lambda x: 1 / x))
    ax2.set_xlabel('Wavelength (km)')
    ax2.set_xticks([2000, 1000, 500, 200, 100, 50, 30, 20])
    ax2.set_xticklabels(['2000', '1000', '500', '200', '100', '50', '30', '20'])

    if ratio:
        ax1.set_xlabel('')
        axr.axhline(1.0, color='0.25', lw=1.1)
        for label, (mean, ci, _, _), color, lw, ls in curves:
            axr.semilogx(master_frequencies, mean / mean_ds, color=color, lw=lw, ls=ls)
        axr.axvspan(1 / MODEL_RES, F_NYQ, color='0.85', alpha=0.5, zorder=0)
        axr.set_yscale('log'); axr.set_ylim(0.04, 4)
        axr.set_yticks([0.05, 0.1, 0.25, 0.5, 1, 2, 4])
        axr.set_yticklabels(['0.05', '0.1', '0.25', '0.5', '1', '2', '4'])
        axr.set_ylabel('model PSD / observed')
        axr.set_xlabel('Wavenumber (cycles / km)')
        axr.grid(True, which='both', ls='--', alpha=0.4)
        axr.text(0.033, 0.06, 'below model\neffective res.', fontsize=8, color='0.35',
                 ha='center', transform=axr.get_xaxis_transform())

    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    print('\nwrote %s' % outfile)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('outfile', nargs='?', default='swot_spectra.png')
    ap.add_argument('-n', '--n-bins', type=int, default=N_BINS)
    ap.add_argument('-m', '--method', choices=['logbin', 'interp'], default=METHOD,
                    help="'interp' reproduces the original np.interp behaviour")
    ap.add_argument('-w', '--weights', choices=['v1', 'v2'], default='v2')
    ap.add_argument('-i', '--iters', nargs=2, type=int, default=[5, 5],
                    metavar=('IT_MG', 'IT_HR'))
    ap.add_argument('-s', '--spec', default=None,
                    help="explicit curve list, e.g. 'MG:v1:8,MG:v2:5,HR:v1:17,HR:v2:10'")
    ap.add_argument('-r', '--ratio', action='store_true',
                    help='add a model/observed ratio panel below')
    a = ap.parse_args()
    mdl = models_from_spec(a.spec) if a.spec else models_for(a.weights, *a.iters)
    main(a.outfile, a.n_bins, a.method, mdl, a.ratio)
