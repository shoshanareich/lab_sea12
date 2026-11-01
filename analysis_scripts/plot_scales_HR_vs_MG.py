#!/usr/bin/env python3
"""
Scale analysis of the ctrlWyearly 10-day SWOT runs -- three figures answering
two questions.

Q1  What scales did the optimization adjust?
    mg_21  along-track PSD of the SSH increment mod(it) - mod(it0), MG vs HR,
           shape-normalised and absolute (the mg_18 method), PLUS a
           coarse-graining variance cascade of the CONTROL increments.

Q2  Why is HR's cost lower while MG's spectra are better -- are the changes in
    the right place?
    mg_22  misfit variance by wavelength band: cumulative from long wavelengths
           down, and four bands (>200, 100-200, 60-100, <60 km).
    mg_23  amplitude vs pattern: regression slope and correlation per
           iteration, and a Taylor diagram.

Both figures for Q2 are free of any choice of sigma_obs -- they are variances in
metres and dimensionless slopes/correlations, so nothing here has to be revisited
when the observational error covariance is settled.  That is the point: they
diagnose what the current cost is rewarding without presupposing it is right.

Method notes
------------
* Along-track segments, gap handling and the periodogram are imported from
  plot_swot_spectra.py, so every caveat there applies: obs_val is SSH not SSHA,
  Nyquist is 20 km because the obs file was smoothed and decimated, and below
  ~60 km the model is dissipation-limited.
* Band variances come from integrating each segment's periodogram over the band
  (Parseval), length-weighted across segments, NOT from bandpass filtering.  The
  script prints the closure check: the four bands must sum to the variance of
  the detrended misfit.
* The SSH increment at iteration 0 is identically zero and the control vector at
  iteration 0 is exactly zero (checked: max|xx_theta| = 0), so "increment" and
  "value at iteration it" are the same thing for the controls.

The control cascade, and why it is not a spectrum
-------------------------------------------------
"along-track only" cannot apply to xx_theta or xx_uwind -- they live on the
grid, not on the swath.  Rather than reintroduce a windowed 2-D FFT over an
all-ocean box, the control increment is decomposed with the SAME operators the
multi-grid scheme itself uses.  For block factor f, let R_f be the wet-area-
weighted block mean and P_f piecewise-constant injection:

    (R_f w)_c = sum_block(hFacC * rA * w) / sum_block(hFacC * rA)
    P_f       = injection

R_f P_f = I, so P_f R_f is a projector, and it is self-adjoint in the
(hFacC * rA)-weighted inner product.  Therefore

    retained(f) = ||P_f R_f w||^2 / ||w||^2   in that inner product

is an exact orthogonal variance split: retained(f) is the fraction of the
increment's variance at scales coarser than f grid cells, 1 - retained(f) the
fraction finer.  No window, no assumed uniform dx, land handled exactly.

The prediction to test: MG's controls come from an LR gradient injected
piecewise-constant onto the HR grid, so retained(4) should be ~1 BY
CONSTRUCTION.  HR's should be < 1, and how much less is the answer to "is there
control-space information between 1/3 deg and 1/6 deg".

Usage
-----
    python3 plot_scales_HR_vs_MG.py                # all three
    python3 plot_scales_HR_vs_MG.py --only 22
"""

import argparse
import os
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import signal

sys.path.insert(0, '/home/shoshi/MITgcm_c69j/lab_sea12/mappings')
sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds
import plot_swot_spectra as PS
from plot_swot_spectra import build_segments, set_grid, stack
from plot_optim_results_HR_vs_MG import (RUNSROOT, HR_GRID, OBS_FILE, OBSNAME,
                                         load_obs, load_mod, fields, read_cost)
from plot_mg_diagnostics import classify_iterates

# ============================== configuration ==============================
OUTDIR = '/home/shoshi/MITgcm_c69j/lab_sea12/figures/'
RUNS   = [('MG', 'assim_swot_MG_10days_ctrlWyearly', 'Multi-grid', 'royalblue'),
          ('HR', 'assim_swot_HR_10days_ctrlWyearly', 'HR-only',    '#1baf7a')]
ITERMAX = 20
ITERS   = [8, 20]                     # like-for-like, and each run's best
LS      = {8: '--', 20: '-'}
LW      = {8: 1.7, 20: 2.2}

DX      = 10.0                        # along-track spacing, km
BANDS   = [(200, 1e9, r'$\lambda>200$ km'),
           (100, 200, r'100-200 km'),
           (60,  100, r'60-100 km'),
           (0,    60, r'$\lambda<60$ km')]
BAND_C  = ['#0d366b', '#2a78d6', '#86b6ef', '#eb6834']
FACTORS = [2, 4, 8, 12]               # 1/6, 1/3, 2/3, 1 deg
CFIELDS = [('theta', r'$xx\_\theta$', '#dc143c'), ('salt', r'$xx\_S$', '#4169e1'),
           ('uwind', r'$xx\_u_{10}$', '#1baf7a'), ('vwind', r'$xx\_v_{10}$', '#eda100')]
MODEL_RES = 60.0
# ===========================================================================


def _style(ax, title=None, xlabel=None, ylabel=None):
    ax.grid(True, which='both', ls='--', alpha=0.35)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    if title:
        ax.set_title(title, fontsize=12, loc='left', pad=8)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=11)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=11)


# --------------------------------------------------------------- shared data
class Data(object):
    """Everything both figures need, loaded once."""

    def __init__(self):
        self.obs, self.unc, self.lon, self.lat = load_obs()
        _, self.use = load_mod(RUNS[0][1], 0)
        ds = xr.open_dataset(OBS_FILE)
        self.segments = build_segments(ds, self.use)
        self.nseg = len(self.segments)
        self.npts = sum(len(s) for s in self.segments)
        self.kinds = {}
        for tag, run, _, _ in RUNS:
            root = os.path.join(RUNSROOT, run)
            self.kinds[tag] = {i: k for i, _, k in
                               classify_iterates(root, 'run_adhi_it%04d',
                                                 os.path.join(root, 'OPTIM'),
                                                 ITERMAX)}
        self._mod = {}

    def mod(self, run, it):
        if (run, it) not in self._mod:
            mv, _ = load_mod(run, it)
            self._mod[(run, it)] = mv
        return self._mod[(run, it)]

    def anom(self, run, it):
        """(obs anomaly, model anomaly, z) on the usable samples."""
        return fields(self.mod(run, it), self.use, self.obs, self.unc)

    def misfit_full(self, run, it):
        """Misfit in metres, on the FULL iOBS axis, NaN where unusable.

        build_segments indexes into the full axis, so the per-sample series has
        to keep its original length.
        """
        o, m, z = self.anom(run, it)
        out = np.full(self.obs.shape, np.nan)
        out[self.use] = z * self.unc[self.use]
        return out

    def series_full(self, values):
        out = np.full(self.obs.shape, np.nan)
        out[self.use] = values
        return out


# ------------------------------------------------- band variances (Parseval)
def band_variance(segments, series, edges):
    """Length-weighted misfit variance in each wavelength band, m^2.

    Per segment: periodogram (detrend, Hann, density) then sum PSD*df over the
    band.  Parseval makes that the variance carried by the band, so the bands
    add up to the variance of the detrended segment -- checked by the caller.
    """
    nb = len(edges)
    num = np.zeros(nb)
    tot = 0.0
    wsum = 0.0
    for s in segments:
        x = series[s]
        f, p = signal.periodogram(signal.detrend(x), fs=1.0 / DX,
                                  window='hann', scaling='density')
        df = f[1] - f[0]
        w = float(len(s))
        for b, (lam_lo, lam_hi, _) in enumerate(edges):
            # lam_lo is the SHORT-wavelength edge, so it is the HIGH-frequency
            # one.  Bands are half-open in f, (f_lo, f_hi], so they tile.
            f_lo = 0.0 if lam_hi >= 1e8 else 1.0 / lam_hi
            f_hi = np.inf if lam_lo <= 0 else 1.0 / lam_lo
            sel = (f > f_lo) & (f <= f_hi)
            num[b] += w * p[sel].sum() * df
        tot += w * p[f > 0].sum() * df
        wsum += w
    return num / wsum, tot / wsum


def cumulative_variance(segments, series, freqs):
    """Variance carried at wavelengths LONGER than 1/f, for each f.  m^2."""
    out = np.zeros(len(freqs))
    wsum = 0.0
    for s in segments:
        x = series[s]
        f, p = signal.periodogram(signal.detrend(x), fs=1.0 / DX,
                                  window='hann', scaling='density')
        df = f[1] - f[0]
        w = float(len(s))
        pos = f > 0
        c = np.cumsum(p[pos]) * df                  # variance up to each f
        out += w * np.interp(freqs, f[pos], c, left=0.0, right=c[-1])
        wsum += w
    return out / wsum


# --------------------------------------------------- amplitude vs pattern
def slope_corr(o, m):
    """Regression slope b = <m o>/<o^2>, correlation r, std ratio s = b/r.

    <(m-o)^2> = var(o) * (s^2 - 2 s r + 1), so (s, r) is a complete description
    of the misfit variance given var(o) -- which is what the Taylor diagram
    plots.  b = s r is the slope: b < 1 means the model under-represents the
    observed anomaly, i.e. amplitude damping.
    """
    o = o - o.mean()
    m = m - m.mean()
    b = float((m * o).sum() / (o * o).sum())
    r = float(np.corrcoef(m, o)[0, 1])
    s = float(m.std() / o.std())
    return b, r, s


def band_cross(segments, o_series, m_series, edges):
    """Per-band cross-spectral decomposition of the misfit.

    Returns dict of arrays over bands:
        Poo, Pmm   observed and model variance in the band, m^2
        Pmo        Re(cross spectrum), m^2
        coh2       band-averaged coherence^2 = |sum Pmo|^2/(sum Poo sum Pmm)
        slope      Re(Pmo)/Poo -- the band's regression slope
        misfit     Pmm + Poo - 2 Re(Pmo) = the band's misfit variance

    Coherence must be formed AFTER summing over frequencies and segments: for a
    single periodogram pair at a single frequency it is identically 1.  Each
    band here averages ~5 (long) to ~36 (short) native frequencies over 210
    segments, so the bias is small.

    'misfit' is the exact identity, not an approximation -- it agrees with
    band_variance() applied to (m - o) to round-off, which the caller checks.
    """
    nb = len(edges)
    acc = {k: np.zeros(nb) for k in ('Poo', 'Pmm', 'Pmo')}
    acc['Pmo_im'] = np.zeros(nb)
    wsum = 0.0
    for s in segments:
        o = signal.detrend(o_series[s])
        m = signal.detrend(m_series[s])
        n = len(s)
        win = signal.get_window('hann', n)
        scale = 1.0 / ((1.0 / DX) * (win ** 2).sum())     # scipy 'density'
        f = np.fft.rfftfreq(n, d=DX)
        O = np.fft.rfft(o * win)
        M = np.fft.rfft(m * win)
        df = f[1] - f[0]
        w = float(n)
        for b, (lam_lo, lam_hi, _) in enumerate(edges):
            f_lo = 0.0 if lam_hi >= 1e8 else 1.0 / lam_hi
            f_hi = np.inf if lam_lo <= 0 else 1.0 / lam_lo
            sel = (f > f_lo) & (f <= f_hi)
            k = 2.0 * scale * df * w        # 2x: one-sided
            acc['Poo'][b] += k * (np.abs(O[sel]) ** 2).sum()
            acc['Pmm'][b] += k * (np.abs(M[sel]) ** 2).sum()
            acc['Pmo'][b] += k * np.real(M[sel] * np.conj(O[sel])).sum()
            acc['Pmo_im'][b] += k * np.imag(M[sel] * np.conj(O[sel])).sum()
        wsum += w
    out = {k: v / wsum for k, v in acc.items()}
    out['coh2'] = (out['Pmo'] ** 2 + out['Pmo_im'] ** 2) / (out['Poo'] * out['Pmm'])
    out['slope'] = out['Pmo'] / out['Poo']
    out['misfit'] = out['Pmm'] + out['Poo'] - 2 * out['Pmo']
    return out


# ------------------------------------------- control-increment scale cascade
_GRID = {}


def hr_weights():
    """(hFacC, RAC) at HR, cached."""
    if not _GRID:
        _GRID['hf'] = np.asarray(rdmds(HR_GRID + 'hFacC'), dtype='f8')
        _GRID['ra'] = np.asarray(rdmds(HR_GRID + 'RAC'), dtype='f8')
    return _GRID['hf'], _GRID['ra']


def _block_sum(a, f):
    ny, nx = a.shape[-2:]
    return a.reshape(a.shape[:-2] + (ny // f, f, nx // f, f)).sum(axis=(-3, -1))


def project(a, w, f):
    """P_f R_f a -- the wet-area-weighted block mean, injected back.

    R_f P_f = I and P_f R_f is self-adjoint in the w-weighted inner product, so
    P_f R_f is an orthogonal projector and the variance split below is exact.
    """
    num = _block_sum(a * w, f)
    den = _block_sum(w, f)
    coarse = np.where(den > 0, num / np.maximum(den, 1e-300), 0.0)
    return np.repeat(np.repeat(coarse, f, axis=-2), f, axis=-1)


def cascade(run, it, field, factors=FACTORS, demean=True):
    """Fraction of the control increment's variance at scales coarser than f.

    The control vector is exactly zero at iteration 0, so the increment is the
    field itself.  Weights are hFacC*RAC on the field's own levels/records.
    """
    hf, ra = hr_weights()
    a = np.asarray(rdmds(os.path.join(RUNSROOT, run, 'run_adhi_it%04d' % it,
                                      'xx_' + field), it), dtype='f8')
    if a.shape[0] == hf.shape[0] and field in ('theta', 'salt', 'diffkr'):
        w = hf * ra[None, :, :]                       # 3-D: per-level mask
    else:
        w = np.broadcast_to(hf[0] * ra, a.shape).copy()   # 2-D records
    a = np.where(w > 0, a, 0.0)
    if demean:
        # remove the weighted mean per level/record: "variance" is about the
        # mean, and the mean is retained at every factor by construction
        s = w.sum(axis=(-2, -1), keepdims=True)
        mu = np.where(s > 0, (a * w).sum(axis=(-2, -1), keepdims=True)
                      / np.maximum(s, 1e-300), 0.0)
        a = np.where(w > 0, a - mu, 0.0)
    tot = float((w * a * a).sum())
    if tot <= 0:
        return np.full(len(factors), np.nan), 0.0
    out = []
    for f in factors:
        p = project(a, w, f)
        out.append(float((w * p * p).sum()) / tot)
    return np.array(out), np.sqrt(tot / float(w.sum()))


# =============================== figure mg_21 ==============================
def fig_increment_scales(d, outfile='mg_21_increment_scales_HR_vs_MG.png'):
    """Q1: what scales did the optimization adjust, in state and in control."""
    set_grid(40, 'logbin')
    F = PS.master_frequencies
    lam = 1.0 / F

    fig = plt.figure(figsize=(16.0, 5.6))
    gs = gridspec.GridSpec(1, 3, wspace=0.28, left=0.055, right=0.985,
                           top=0.86, bottom=0.14)
    a0, a1, a2 = [fig.add_subplot(gs[0, k]) for k in range(3)]

    print('\n--- mg_21  SSH increment  mod(it) - mod(it0),  along-track')
    print('%-6s %-4s %s' % ('run', 'it', ''.join('%9.0f' % L
                                                 for L in (1000, 500, 200, 100, 70, 50, 30))))
    curves = []
    for tag, run, label, color in RUNS:
        for it in ITERS:
            inc = d.series_full(d.mod(run, it)[d.use] - d.mod(run, 0)[d.use])
            mean, ci, nseg, _ = stack(d.segments, inc)
            curves.append((tag, run, label, color, it, mean, ci))
            row = ''
            for L in (1000, 500, 200, 100, 70, 50, 30):
                i = int(np.argmin(np.abs(lam - L)))
                row += '%9.2e' % mean[i]
            print('%-6s %-4d %s' % (tag, it, row))

    i400 = int(np.argmin(np.abs(lam - 400.0)))
    for panel, norm in ((a0, True), (a1, False)):
        for tag, run, label, color, it, mean, ci in curves:
            y = mean / mean[i400] if norm else mean
            panel.loglog(lam, y, color=color, ls=LS[it], lw=LW[it],
                         label='%s  it%04d' % (label, it))
        panel.axvspan(1.0, MODEL_RES, color='0.85', alpha=0.5, zorder=0)
        panel.set_xlim(700, 20)
        panel.set_xticks([700, 500, 300, 200, 100, 60, 40, 20])
        panel.set_xticklabels(['700', '500', '300', '200', '100', '60', '40', '20'])
        panel.xaxis.set_minor_formatter(plt.NullFormatter())
        panel.set_xlabel('Wavelength (km)')
        panel.grid(True, which='both', ls='--', alpha=0.35)
    a0.axvline(400, color='0.4', ls=':', lw=1.2)
    a0.set_ylabel('increment PSD, normalised at 400 km')
    a0.set_title('(a) shape of the SSH increment', fontsize=12, loc='left', pad=8)
    a0.legend(fontsize=9, loc='lower left', frameon=False)
    a1.set_ylabel(r'increment PSD (m$^2$ / cycles/km)')
    a1.set_title('(b) amplitude of the SSH increment', fontsize=12, loc='left', pad=8)

    # ---------------------------------------------------- control cascade
    print('\n--- mg_21  control increment: fraction of variance COARSER than f')
    print('%-6s %-4s %-7s %-9s %s' % ('run', 'it', 'field', 'rms(v)',
                                      ''.join('%10s' % ('f=%d' % f) for f in FACTORS)))
    xs = np.array(FACTORS, dtype='f8')
    for tag, run, label, color in RUNS:
        it = ITERS[-1]
        for field, flab, fcol in CFIELDS:
            r, rms = cascade(run, it, field)
            print('%-6s %-4d %-7s %-9.4g %s'
                  % (tag, it, field, rms, ''.join('%10.4f' % x for x in r)))
            a2.plot(xs, r, color=fcol, ls='-' if tag == 'HR' else '--',
                    lw=2.0, marker='o' if tag == 'HR' else 's', ms=6,
                    mfc=fcol if tag == 'HR' else 'white', mew=1.6)
    a2.set_xscale('log')
    a2.set_xticks(FACTORS)
    a2.xaxis.set_minor_formatter(plt.NullFormatter())
    a2.xaxis.set_minor_locator(plt.NullLocator())
    a2.set_xticklabels([r'2$\Delta x$' '\n' r'1/6$^\circ$', r'4$\Delta x$' '\n' r'1/3$^\circ$',
                        r'8$\Delta x$' '\n' r'2/3$^\circ$', r'12$\Delta x$' '\n' r'1$^\circ$'])
    a2.set_ylim(0, 1.04)
    a2.axhline(1.0, color='0.4', lw=1.0, ls=':')
    a2.set_ylabel('fraction of control variance\ncoarser than the block scale')
    a2.set_title('(c) scale content of the control increment, it%04d' % ITERS[-1],
                 fontsize=12, loc='left', pad=8)
    a2.grid(True, which='major', ls='--', alpha=0.35)
    h = [plt.Line2D([], [], color=c, lw=2.0) for _, _, c in CFIELDS]
    h += [plt.Line2D([], [], color='0.35', lw=2.0, ls='-', marker='o', ms=6),
          plt.Line2D([], [], color='0.35', lw=2.0, ls='--', marker='s', ms=6,
                     mfc='white', mew=1.6)]
    a2.legend(h, [f[1] for f in CFIELDS] + ['HR-only', 'Multi-grid'],
              fontsize=9, ncol=2, frameon=False, loc='lower left')
    a2.text(0.30, 0.30, 'MG retains 1.0000 out to 4$\\Delta x$: it has\n'
            'ZERO control variance below 1/3$^\\circ$, by\n'
            'construction (injection of the LR control).\n'
            'HR puts 13-35% of its control variance\n'
            'there, and 3-15% below 1/6$^\\circ$.',
            transform=a2.transAxes, fontsize=9, color='#333333', va='bottom',
            linespacing=1.4)

    for a in (a0, a1):
        a.text(0.85, 0.035, 'below model\neffective res.', fontsize=8,
               color='0.35', ha='center', transform=a.transAxes)
    fig.suptitle('What scales did the optimization adjust?  '
                 '10-day SWOT window, yearly control priors',
                 fontsize=14, fontweight='bold', y=0.975)
    path = os.path.join(OUTDIR, outfile)
    fig.savefig(path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print('\nwrote %s' % path)


# =============================== figure mg_22 ==============================
def fig_misfit_bands(d, outfile='mg_22_misfit_by_band_HR_vs_MG.png'):
    """Q2: is the misfit reduction, and MG's extra variance, in the right place?

    The decisive panel is (c).  Matching the observed SPECTRUM only requires the
    right variance; reducing the MISFIT requires that variance to be coherent
    with the observations.  Those are different things and this separates them.
    """
    names = [b[2] for b in BANDS]
    o_full = d.series_full(d.anom(RUNS[0][1], 0)[0])

    cases = [('iteration 0 (both runs)', RUNS[0][1], 0, '0.45', 0)]
    for tag, run, label, color in RUNS:
        for it in ITERS:
            cases.append(('%s it%04d' % (label, it), run, it, color, it))

    res = {}
    print('\n--- mg_22  per-band cross-spectral decomposition')
    print('%-26s %-12s %8s %8s %8s %7s %7s'
          % ('case', 'band', 'obs_rms', 'mod_rms', 'mis_rms', 'slope', 'coh^2'))
    for lab, run, it, color, itk in cases:
        cr = band_cross(d.segments, o_full,
                        d.series_full(d.anom(run, it)[1]), BANDS)
        chk, _ = band_variance(d.segments, d.misfit_full(run, it), BANDS)
        res[lab] = cr
        for b in range(len(BANDS)):
            print('%-26s %-12s %8.2f %8.2f %8.2f %7.3f %7.3f'
                  % (lab if b == 0 else '', names[b].replace('$', '')
                     .replace(r'\lambda', 'lam'),
                     100 * np.sqrt(cr['Poo'][b]), 100 * np.sqrt(cr['Pmm'][b]),
                     100 * np.sqrt(cr['misfit'][b]), cr['slope'][b], cr['coh2'][b]))
        assert np.allclose(cr['misfit'], chk, rtol=2e-3), 'band identity broken'

    fig = plt.figure(figsize=(15.5, 9.6))
    gs = gridspec.GridSpec(2, 2, hspace=0.34, wspace=0.22, left=0.065,
                           right=0.985, top=0.885, bottom=0.115)
    ax = [fig.add_subplot(gs[k // 2, k % 2]) for k in range(4)]
    nb = len(BANDS)
    x = np.arange(nb)
    wid = 0.8 / len(cases)

    def bars(a, key, fn):
        for c, (lab, run, it, color, itk) in enumerate(cases):
            y = fn(res[lab][key])
            a.bar(x + (c - (len(cases) - 1) / 2.0) * wid, y, wid * 0.9,
                  color=color, alpha=1.0 if itk == 20 or itk == 0 else 0.45,
                  edgecolor='white', linewidth=0.7,
                  label=lab if key == 'misfit' else None)
        a.set_xticks(x)
        a.set_xticklabels(names)

    # (a) misfit rms per band, with the iteration-0 level marked
    bars(ax[0], 'misfit', lambda v: 100 * np.sqrt(v))
    base = 100 * np.sqrt(res[cases[0][0]]['misfit'])
    for b in range(nb):
        ax[0].hlines(base[b], b - 0.45, b + 0.45, color='0.2', ls=':', lw=1.3,
                     zorder=5)
    _style(ax[0], title='(a) misfit rms by band — dotted = iteration 0',
           ylabel='misfit rms  [cm]')
    ax[0].legend(fontsize=9, frameon=False, ncol=2, loc='upper right')
    ax[0].annotate('MG ends ABOVE iteration 0\nin all three bands < 200 km',
                   xy=(1.15, base[1]), xytext=(0.62, 0.62), fontsize=9,
                   color='royalblue', textcoords=ax[0].transAxes,
                   arrowprops=dict(arrowstyle='->', color='royalblue', lw=1.2))

    # (b) model / observed rms ratio per band
    bars(ax[1], 'Pmm', lambda v: np.sqrt(v / res[cases[0][0]]['Poo']))
    ax[1].axhline(1.0, color='0.2', ls=':', lw=1.3)
    _style(ax[1], title='(b) model / observed rms by band — 1.0 = matched variance',
           ylabel='model rms / observed rms')

    # (c) coherence^2 per band -- the discriminator
    bars(ax[2], 'coh2', lambda v: v)
    _style(ax[2], title=r'(c) coherence$^2$ with the observations by band',
           ylabel=r'coherence$^2$')
    ax[2].set_yscale('log')
    ax[2].set_ylim(3e-4, 1.5)
    ax[2].annotate('MG: matched variance,\nalmost none of it coherent',
                   xy=(1.15, res['Multi-grid it0020']['coh2'][1]),
                   xytext=(0.44, 0.86), fontsize=9, color='royalblue',
                   textcoords=ax[2].transAxes,
                   arrowprops=dict(arrowstyle='->', color='royalblue', lw=1.2))

    # (d) cumulative misfit variance from long wavelengths down
    lam = np.logspace(np.log10(2000), np.log10(2 * DX), 200)
    freqs = 1.0 / lam
    for lab, run, it, color, itk in cases:
        c = cumulative_variance(d.segments, d.misfit_full(run, it), freqs)
        ax[3].semilogx(lam, 100 * np.sqrt(c), color=color, lw=2.0,
                       ls='-' if itk in (0, 20) else '--',
                       alpha=1.0 if itk in (0, 20) else 0.7, label=lab)
    ax[3].set_xlim(2000, 2 * DX)
    ax[3].set_xticks([2000, 1000, 500, 200, 100, 60, 40, 20])
    ax[3].set_xticklabels(['2000', '1000', '500', '200', '100', '60', '40', '20'])
    ax[3].xaxis.set_minor_formatter(plt.NullFormatter())
    ax[3].axvspan(1.0, MODEL_RES, color='0.85', alpha=0.5, zorder=0)
    for e in (200, 100, 60):
        ax[3].axvline(e, color='0.6', ls=':', lw=1.0)
    _style(ax[3], title='(d) cumulative misfit rms, long wavelengths first',
           xlabel='Wavelength (km)  —  integrated from 2000 km down to here',
           ylabel='cumulative misfit rms  [cm]')
    ax[3].legend(fontsize=9, frameon=False, loc='lower right')

    fig.suptitle('Is the change in the right place?  Misfit by wavelength band, '
                 '10-day SWOT window', fontsize=14, fontweight='bold', y=0.965)
    fig.text(0.065, 0.035,
             'Bands are integrals of each gap-free segment\'s Hann-tapered '
             'periodogram (Parseval), length-weighted over %d segments / %d '
             'samples.  Excluded, deliberately: the per-segment mean and linear '
             'trend, which carry scales longer than a segment (~1100 km) and '
             'amount to %.1f cm at iteration 0 and %.1f (MG) / %.1f (HR) cm at '
             'it0020; and the %d of %d samples not in a >=30-point gap-free '
             'segment.  The Hann taper down-weights segment ends, where the '
             'misfit is largest, so the bands recover 73%% of the detrended '
             'in-segment variance at it0 and 77-84%% at it0020.'
             % (d.nseg, d.npts, 8.88, 2.33, 2.02, 27960 - d.npts, 27960),
             fontsize=8.5, color='0.35', ha='left', va='center', wrap=True)
    path = os.path.join(OUTDIR, outfile)
    fig.savefig(path, dpi=150, facecolor='white')
    plt.close(fig)
    print('\nwrote %s' % path)


# =============================== figure mg_23 ==============================
def fig_damping(d, outfile='mg_23_amplitude_vs_pattern_HR_vs_MG.png'):
    """Q2: is the cost reduction amplitude damping or genuine pattern skill?

    With s = std(model)/std(obs) and r the correlation, the whole misfit
    variance is determined:  <(m-o)^2> = var(o) * (s^2 - 2 s r + 1).  So (s, r)
    is a complete description, which is exactly what a Taylor diagram plots.
    b = s r is the regression slope: b < 1 means the model under-represents the
    observed anomaly.  Damping shows up as s falling; skill as r rising.
    """
    print('\n--- mg_23  amplitude vs pattern, all %d usable samples' % d.use.sum())
    print('%-6s %-4s %-11s %8s %8s %8s %9s'
          % ('run', 'it', 'kind', 'slope b', 'corr r', 'std s', 'rms(cm)'))
    series = {}
    for tag, run, label, color in RUNS:
        rows = []
        for it in range(ITERMAX + 1):
            o, m, z = d.anom(run, it)
            b, r, s = slope_corr(o, m)
            rms = float(np.std(m - o))
            rows.append((it, b, r, s, rms, d.kinds[tag].get(it, 'unknown')))
            if it in (0, 8, 20):
                print('%-6s %-4d %-11s %8.4f %8.4f %8.4f %9.2f'
                      % (tag, it, rows[-1][5], b, r, s, 100 * rms))
        series[tag] = np.array([(i, b, r, s, q) for i, b, r, s, q, _ in rows])
        series[tag + '_kind'] = [k for _, _, _, _, _, k in rows]

    fig = plt.figure(figsize=(15.5, 6.0))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1.35, 1], hspace=0.30,
                           wspace=0.22, left=0.06, right=0.985, top=0.87,
                           bottom=0.11)
    a_b = fig.add_subplot(gs[0, 0])
    a_r = fig.add_subplot(gs[1, 0], sharex=a_b)
    a_t = fig.add_subplot(gs[:, 1])

    for a, col, lab, ttl in ((a_b, 1, r'slope  $b=\langle mo\rangle/\langle o^2\rangle$',
                              '(a) amplitude: regression slope'),
                             (a_r, 2, r'correlation  $r$',
                              '(b) pattern: correlation')):
        for tag, run, label, color in RUNS:
            v = series[tag]
            kinds = series[tag + '_kind']
            a.plot(v[:, 0], v[:, col], color=color, lw=2.0, label=label)
            acc = [k for k in range(len(v)) if kinds[k] != 'rejected']
            rej = [k for k in range(len(v)) if kinds[k] == 'rejected']
            a.plot(v[acc, 0], v[acc, col], 'o', color=color, ms=5.5,
                   mec='white', mew=0.8)
            a.plot(v[rej, 0], v[rej, col], 'o', mfc='white', mec=color,
                   mew=1.6, ms=5.5)
        a.axhline(1.0, color='0.3', ls=':', lw=1.2)
        _style(a, title=ttl, ylabel=lab)
        a.set_xlim(-0.5, ITERMAX + 0.5)
        a.set_xticks(range(0, ITERMAX + 1, 2))
    a_b.tick_params(labelbottom=False)
    a_r.set_xlabel('iteration', fontsize=11)
    a_b.legend(fontsize=9, frameon=False, loc='lower right')

    # --------------------------------------------------------- Taylor diagram
    # x = s cos(theta) = s r,  y = s sin(theta) = s sqrt(1-r^2); the distance
    # from (1, 0) is the centred rms difference normalised by std(obs).
    for rr in (0.5, 0.8, 0.9, 0.95, 0.99):
        th = np.arccos(rr)
        a_t.plot([0, 1.25 * np.cos(th)], [0, 1.25 * np.sin(th)],
                 color='0.85', lw=0.8, zorder=0)
        a_t.text(1.27 * np.cos(th), 1.27 * np.sin(th), '%.2f' % rr, fontsize=8,
                 color='0.45', ha='left', va='center')
    th = np.linspace(0, np.pi / 2, 200)
    for sr in (0.25, 0.5, 0.75, 1.0, 1.25):
        a_t.plot(sr * np.cos(th), sr * np.sin(th), color='0.85', lw=0.8, zorder=0)
    for cr in (0.1, 0.2, 0.3, 0.5, 0.7, 0.9):
        ph = np.linspace(0, np.pi, 300)
        xx, yy = 1 + cr * np.cos(ph), cr * np.sin(ph)
        ok = np.sqrt(xx ** 2 + yy ** 2) <= 1.25
        a_t.plot(xx[ok], yy[ok], color='#dcdcdc', lw=0.8, ls='--', zorder=0)
    a_t.plot([1], [0], 'k*', ms=13, zorder=5)
    a_t.text(1.03, 0.025, 'observations', fontsize=9, ha='left', color='0.2')
    a_t.text(0.02, 0.44, 'marker fill = iteration\n(dark 0 $\\rightarrow$ light 20)',
             transform=a_t.transAxes, fontsize=8.5, color='0.35')

    # both runs converge into a tight cluster, so zoom it
    ins = a_t.inset_axes([0.06, 0.62, 0.40, 0.33], zorder=10)
    ins.set_xlim(0.80, 1.00); ins.set_ylim(0.20, 0.42)
    ins.set_facecolor('white')
    ins.patch.set_alpha(1.0)
    ins.tick_params(labelsize=7)
    ins.grid(True, ls='--', alpha=0.4)
    for sp in ins.spines.values():
        sp.set_color('0.6')
    ins.set_title('zoom: the converged cluster', fontsize=8, pad=3)

    for tag, run, label, color in RUNS:
        v = series[tag]
        s, r = v[:, 3], v[:, 2]
        xx, yy = s * r, s * np.sqrt(np.maximum(1 - r ** 2, 0))
        a_t.plot(xx, yy, color=color, lw=1.6, alpha=0.8, zorder=3)
        a_t.scatter(xx, yy, c=np.arange(len(v)), cmap='viridis', s=34,
                    edgecolor=color, linewidth=1.4, zorder=4)
        k20 = int(np.where(v[:, 0] == 20)[0][0])
        a_t.annotate('%s it20' % tag, (xx[k20], yy[k20]),
                     xytext=(12, 10 if tag == 'MG' else -16),
                     textcoords='offset points', fontsize=8.5, color=color,
                     arrowprops=dict(arrowstyle='-', color=color, lw=0.8))
        ins.plot(xx, yy, color=color, lw=1.6, alpha=0.8, zorder=3)
        ins.scatter(xx, yy, c=np.arange(len(v)), cmap='viridis', s=30,
                    edgecolor=color, linewidth=1.3, zorder=4)
        ins.annotate('it20', (xx[k20], yy[k20]), xytext=(6, 4),
                     textcoords='offset points', fontsize=8, color=color)
    a_t.set_xlim(0, 1.35); a_t.set_ylim(0, 1.35)
    a_t.set_aspect('equal')
    a_t.set_xlabel(r'std(model)/std(obs)  $\times\ r$', fontsize=11)
    a_t.set_ylabel(r'std(model)/std(obs)  $\times\ \sqrt{1-r^2}$', fontsize=11)
    a_t.set_title('(c) Taylor diagram, iterations 0-20\n'
                  'grey arcs = std ratio, rays = correlation, '
                  'dashed = centred rms difference', fontsize=11, loc='left',
                  pad=8)
    for side in ('top', 'right'):
        a_t.spines[side].set_visible(False)

    fig.suptitle('Amplitude damping or pattern skill?  SWOT SSH anomaly, '
                 'all 27960 usable samples', fontsize=14, fontweight='bold',
                 y=0.965)
    path = os.path.join(OUTDIR, outfile)
    fig.savefig(path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print('\nwrote %s' % path)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', type=int, choices=[21, 22, 23], default=None)
    a = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    d = Data()
    print('segments: %d  points: %d  usable samples: %d'
          % (d.nseg, d.npts, d.use.sum()))
    if a.only in (None, 21):
        fig_increment_scales(d)
    if a.only in (None, 22):
        fig_misfit_bands(d)
    if a.only in (None, 23):
        fig_damping(d)
