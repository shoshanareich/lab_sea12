#!/usr/bin/env python3
"""
Three plain-language figures for one conclusion:

    the multi-grid run's better SSH spectrum is variance without coherence.

mg_24  Correlation with SWOT per outer iteration.  MG saturates at ~0.95 by
       iteration 4 and never improves; HR keeps climbing to 0.969.  Second
       panel is the consequence in cm.
mg_25  ONE along-track segment.  Top: the SSH along it.  Bottom: the same
       segment band-passed to 60-200 km, where MG carries 98% of the observed
       amplitude and HR only 53% -- yet HR is the one that lines up.
mg_26  Two bar charts.  Left: model/observed variance by band, where MG sits at
       1.0 and looks ideal.  Right: the fraction of observed variance the model
       actually EXPLAINS (coherence^2), where MG is ~0.  Same runs, opposite
       verdicts, which is the whole point.

Caveat carried in mg_25's caption: a single 167-point segment holds only ~5-15
independent Fourier modes in the 60-200 km band, so its own correlation is
biased high (0.41 MG / 0.64 HR here against population coherence of 0.22 / 0.42).
The segment is an ILLUSTRATION; the number to quote is the population coherence
in mg_26, which averages 210 segments.

Everything here is in cm or dimensionless, so none of it depends on the choice
of sigma_obs.

Usage:  python3 plot_takehome_HR_vs_MG.py [--only 24|25|26]
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, '/home/shoshi/MITgcm_c69j/lab_sea12/mappings')
import plot_scales_HR_vs_MG as S
from plot_scales_HR_vs_MG import (Data, RUNS, ITERMAX, BANDS, OUTDIR, DX,
                                  band_cross, slope_corr, _style)

SEG   = 50          # 167 points; MG amplitude ratio 0.98, HR 0.53 in-band
BAND  = (60.0, 200.0)
OBSC  = 'black'
# ===========================================================================


def bandpass(x, lo=BAND[0], hi=BAND[1], dx=DX):
    """Zero out everything outside [1/hi, 1/lo].  Trim the ends when plotting."""
    x = x - x.mean()
    n = len(x)
    f = np.fft.rfftfreq(n, d=dx)
    X = np.fft.rfft(x)
    return np.fft.irfft(X * ((f >= 1.0 / hi) & (f <= 1.0 / lo)), n=n)


# =============================== figure mg_24 ==============================
def fig_correlation(d, outfile='mg_24_correlation_by_iteration.png'):
    print('\n--- mg_24  correlation and misfit per iteration')
    ser = {}
    for tag, run, label, color in RUNS:
        rows = []
        for it in range(ITERMAX + 1):
            o, m, z = d.anom(run, it)
            b, r, s = slope_corr(o, m)
            rows.append((it, r, 100 * float(np.std(m - o))))
        ser[tag] = np.array(rows)
        ser[tag + 'k'] = [d.kinds[tag].get(i, 'unknown') for i in range(ITERMAX + 1)]
        print('  %-11s r: it0 %.3f  it4 %.3f  it8 %.3f  it20 %.3f   |  '
              'rms cm: it0 %.1f  it20 %.1f'
              % (label, ser[tag][0, 1], ser[tag][4, 1], ser[tag][8, 1],
                 ser[tag][20, 1], ser[tag][0, 2], ser[tag][20, 2]))

    fig, ax = plt.subplots(2, 1, figsize=(9.6, 7.4), sharex=True,
                           gridspec_kw=dict(hspace=0.13, height_ratios=[1.35, 1]))
    for a, col, ylab, ttl in ((ax[0], 1, 'correlation with SWOT SSH',
                               '(a) how well the model matches the OBSERVED PATTERN'),
                              (ax[1], 2, 'misfit rms  [cm]',
                               '(b) and what that is worth, in cm')):
        for tag, run, label, color in RUNS:
            v, kk = ser[tag], ser[tag + 'k']
            a.plot(v[:, 0], v[:, col], color=color, lw=2.4, label=label, zorder=3)
            acc = [i for i in range(len(v)) if kk[i] != 'rejected']
            rej = [i for i in range(len(v)) if kk[i] == 'rejected']
            a.plot(v[acc, 0], v[acc, col], 'o', color=color, ms=6.5, mec='white',
                   mew=1.0, zorder=4)
            a.plot(v[rej, 0], v[rej, col], 'o', mfc='white', mec=color, mew=1.8,
                   ms=6.5, zorder=4)
        _style(a, title=ttl, ylabel=ylab)
        a.set_xlim(-0.5, ITERMAX + 0.5)
        a.set_xticks(range(0, ITERMAX + 1, 2))
    ax[1].set_xlabel('outer iteration', fontsize=11)

    # the point of the figure, said on the figure
    mg, hr = ser['MG'], ser['HR']
    ax[0].axvline(4, color='royalblue', ls=':', lw=1.4, zorder=1)
    ax[0].annotate('MG stops improving here\n(r = %.3f at it4, %.3f at it20 —\n'
                   '16 more iterations, no gain in pattern)'
                   % (mg[4, 1], mg[20, 1]),
                   xy=(4, mg[4, 1]), xytext=(6.4, 0.74), fontsize=9.5,
                   color='royalblue',
                   arrowprops=dict(arrowstyle='->', color='royalblue', lw=1.3))
    ax[0].annotate('HR keeps going: %.3f' % hr[20, 1], xy=(20, hr[20, 1]),
                   xytext=(15.0, 0.855), fontsize=9.5, color='#1baf7a',
                   arrowprops=dict(arrowstyle='->', color='#1baf7a', lw=1.3))
    h, l = ax[0].get_legend_handles_labels()
    h += [plt.Line2D([], [], color='0.35', ls='', marker='o', ms=6.5),
          plt.Line2D([], [], mfc='white', mec='0.35', mew=1.8, ls='', marker='o',
                     ms=6.5)]
    ax[0].legend(h, [r[2] for r in RUNS] + ['accepted iterate', 'rejected trial'],
                 fontsize=9.5, frameon=False, loc='lower right')

    fig.suptitle('Multi-grid stops learning the pattern at iteration 4',
                 fontsize=15, fontweight='bold', y=0.965)
    fig.text(0.5, 0.012, 'Correlation and misfit over all 27960 usable SWOT '
             'samples.  Both runs are bit-identical at iteration 0.',
             fontsize=9, color='0.4', ha='center')
    path = os.path.join(OUTDIR, outfile)
    fig.savefig(path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print('wrote %s' % path)


# =============================== figure mg_25 ==============================
def fig_one_segment(d, seg=SEG, outfile='mg_25_one_segment_bandpassed.png'):
    """The whole conclusion on one satellite track."""
    idx = d.segments[seg]
    o_a = d.anom(RUNS[0][1], 0)[0]
    o = d.series_full(o_a)[idx]
    mg = d.series_full(d.anom('assim_swot_MG_10days_ctrlWyearly', 20)[1])[idx]
    hr = d.series_full(d.anom('assim_swot_HR_10days_ctrlWyearly', 20)[1])[idx]
    dist = np.arange(len(idx)) * DX
    lo, la = d.lon[idx], d.lat[idx]

    t = slice(6, -6)                      # trim FFT edge ringing before plotting
    O, M, H = bandpass(o)[t], bandpass(mg)[t], bandpass(hr)[t]
    D = dist[t]
    print('\n--- mg_25  segment %d: %d points, %.0f km, '
          'lon %.1f-%.1f E, lat %.1f-%.1f N'
          % (seg, len(idx), dist[-1], lo.min(), lo.max(), la.min(), la.max()))
    for nm, X in (('SWOT', O), ('MG', M), ('HR', H)):
        r = np.corrcoef(X, O)[0, 1]
        print('  %-5s band 60-200 km: rms %.2f cm  (%.0f%% of observed)  r = %.2f'
              % (nm, 100 * X.std(), 100 * X.std() / O.std(), r))

    fig, ax = plt.subplots(2, 1, figsize=(12.4, 8.0), sharex=True,
                           gridspec_kw=dict(hspace=0.16, height_ratios=[1, 1.25]))

    ax[0].plot(dist, 100 * (o - o.mean()), color=OBSC, lw=2.6, label='SWOT')
    ax[0].plot(dist, 100 * (mg - mg.mean()), color='royalblue', lw=1.8,
               label='Multi-grid it0020')
    ax[0].plot(dist, 100 * (hr - hr.mean()), color='#1baf7a', lw=1.8,
               label='HR-only it0020')
    _style(ax[0], title='(a) SSH along one SWOT track — at this scale all three agree',
           ylabel='SSH anomaly  [cm]')
    ax[0].legend(fontsize=10, frameon=False, ncol=1, loc='upper left')

    ax[1].axhline(0, color='0.75', lw=1.0)
    ax[1].plot(D, 100 * O, color=OBSC, lw=2.8, zorder=4,
               label='SWOT                       %.1f cm' % (100 * O.std()))
    ax[1].plot(D, 100 * M, color='royalblue', lw=2.2, zorder=3,
               label='Multi-grid    %3.0f%% as big, r = %.2f'
                     % (100 * M.std() / O.std(), np.corrcoef(M, O)[0, 1]))
    ax[1].plot(D, 100 * H, color='#1baf7a', lw=2.2, zorder=3,
               label='HR-only        %3.0f%% as big, r = %.2f'
                     % (100 * H.std() / O.std(), np.corrcoef(H, O)[0, 1]))
    _style(ax[1], title='(b) the same track, keeping only 60-200 km — '
                        'MG is the right SIZE, HR is in the right PLACE',
           xlabel='distance along track  [km]', ylabel='SSH, 60-200 km band  [cm]')
    ax[1].legend(fontsize=10.5, frameon=False, loc='upper right',
                 prop=dict(family='DejaVu Sans Mono', size=9.5))

    fig.suptitle('Matching the amount of wiggle is not the same as matching '
                 'the wiggles', fontsize=15, fontweight='bold', y=0.965)
    fig.text(0.5, 0.008,
             'Segment %d of 210: %d gap-free samples, %.0f km, %.1f-%.1f$^\\circ$N.  '
             'A single segment holds only ~5-15 independent Fourier modes in this '
             'band, so its own r is biased high — the numbers to quote are the '
             '210-segment coherences in mg_26.'
             % (seg, len(idx), dist[-1], la.min(), la.max()),
             fontsize=9, color='0.4', ha='center')
    path = os.path.join(OUTDIR, outfile)
    fig.savefig(path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print('wrote %s' % path)


# =============================== figure mg_26 ==============================
def fig_two_verdicts(d, outfile='mg_26_variance_vs_coherence.png'):
    """Same two runs, two bar charts, opposite conclusions."""
    names = [b[2] for b in BANDS]
    o_full = d.series_full(d.anom(RUNS[0][1], 0)[0])
    cr = {}
    for tag, run, label, color in RUNS:
        cr[tag] = band_cross(d.segments, o_full,
                             d.series_full(d.anom(run, 20)[1]), BANDS)

    fig, ax = plt.subplots(1, 2, figsize=(14.2, 5.8))
    x = np.arange(len(BANDS))
    w = 0.36
    for k, (tag, run, label, color) in enumerate(RUNS):
        off = (k - 0.5) * w
        ax[0].bar(x + off, np.sqrt(cr[tag]['Pmm'] / cr[tag]['Poo']), w * 0.92,
                  color=color, edgecolor='white', linewidth=0.8, label=label)
        ax[1].bar(x + off, cr[tag]['coh2'], w * 0.92, color=color,
                  edgecolor='white', linewidth=0.8, label=label)
    for a in ax:
        a.set_xticks(x)
        a.set_xticklabels(names, fontsize=11)
        a.set_ylim(0, 1.30)
    ax[0].axhline(1.0, color='0.2', ls=':', lw=1.6)
    ax[0].text(3.45, 1.02, 'perfect', fontsize=9, color='0.3', ha='right')
    _style(ax[0], title='(a) Does the model have the right AMOUNT of variance?',
           ylabel='model rms / observed rms')
    _style(ax[1], title='(b) Is that variance in the right PLACE?',
           ylabel=r'fraction of observed variance explained  (coherence$^2$)')
    ax[0].legend(fontsize=10.5, frameon=False, loc='lower left')

    for k, (tag, run, label, color) in enumerate(RUNS):
        off = (k - 0.5) * w
        for b in range(len(BANDS)):
            ax[0].text(b + off, np.sqrt(cr[tag]['Pmm'][b] / cr[tag]['Poo'][b]) + 0.02,
                       '%.2f' % np.sqrt(cr[tag]['Pmm'][b] / cr[tag]['Poo'][b]),
                       ha='center', fontsize=8.5, color=color)
            ax[1].text(b + off, cr[tag]['coh2'][b] + 0.02,
                       '%.2f' % cr[tag]['coh2'][b], ha='center', fontsize=8.5,
                       color=color)

    ax[0].annotate('MG looks IDEAL here: it matches the\n'
                   'observed variance at every scale',
                   xy=(1.18, 1.07), xytext=(0.40, 0.94),
                   textcoords=ax[0].transAxes, fontsize=10, color='royalblue',
                   zorder=6, arrowprops=dict(arrowstyle='->', color='royalblue',
                                             lw=1.4))
    ax[1].annotate('...and explains essentially\nNONE of it below 200 km.\n'
                   'HR explains 3-7x more.',
                   xy=(2, 0.03), xytext=(0.36, 0.44), textcoords=ax[1].transAxes,
                   fontsize=10, color='royalblue', zorder=6,
                   arrowprops=dict(arrowstyle='->', color='royalblue', lw=1.4))

    fig.suptitle('The same two runs, judged two ways — both at their best '
                 'iterate (it0020)', fontsize=15, fontweight='bold', y=0.97)
    fig.text(0.5, 0.005, 'Cross-spectra over 210 gap-free along-track segments '
             '/ 23445 samples.  A spectrum only tests (a).  Only (b) reduces '
             'the misfit.', fontsize=9.5, color='0.4', ha='center')
    fig.tight_layout(rect=[0, 0.035, 1, 0.93])
    path = os.path.join(OUTDIR, outfile)
    fig.savefig(path, dpi=150, facecolor='white')
    plt.close(fig)
    print('wrote %s' % path)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', type=int, choices=[24, 25, 26], default=None)
    a = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    d = Data()
    if a.only in (None, 24):
        fig_correlation(d)
    if a.only in (None, 25):
        fig_one_segment(d)
    if a.only in (None, 26):
        fig_two_verdicts(d)
