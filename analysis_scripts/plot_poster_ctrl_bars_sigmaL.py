#!/usr/bin/env python3
"""
Poster figure: how much each control moved, HR-only vs MG, 10-day sigmaL runs.

    ctrl_adjustment_bars_HR_vs_MG_sigmaL_iter10

RMS over wet points of the raw dimensionless control v (xx_<name> in the HR run
dir, iteration 10).  The penalty is sum(v^2) under WC01, so v is in units of the
prior sigma and the controls are directly comparable; plotted as % of sigma.  For MG these are the
LR controls after prolongation to 1/12 deg, i.e. what the HR forward applied.
3-D controls (theta, salt, diffkr) are RMS over all wet levels; 2-D forcing over
all wet surface points and all daily records.

Usage:  python3 plot_poster_ctrl_bars_sigmaL.py [--png-dir DIR] [--linear]
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from plot_poster_maps_sigmaL import rdmds, native, HR_RUN, MG_RUN, HR_GRID, OUTDIR

ITER = 10
NAME = 'ctrl_adjustment_bars_HR_vs_MG_sigmaL_iter%d' % ITER
CTRLS = [('theta', 'Initial temperature'), ('salt', 'Initial salinity'),
         ('uwind', 'Zonal wind'), ('vwind', 'Meridional wind'),
         ('atemp', 'Air temperature'), ('aqh', 'Humidity'),
         ('precip', 'Precipitation'), ('lwdown', 'Longwave'),
         ('swdown', 'Shortwave'), ('diffkr', 'Vertical mixing')]
C_MG, C_HR = 'royalblue', '#1baf7a'   # as in plot_cost_HR_vs_MG_sigmaL.py
FLOOR = 1e-3     # percent; bars start here on the log axis


def rms_v(run, ctrl, wet3):
    a = native(rdmds(run + 'run_adhi_it%04d/xx_%s' % (ITER, ctrl), itrs=np.nan))
    m = wet3 if a.shape == wet3.shape else np.broadcast_to(wet3[0], a.shape)
    return np.sqrt(np.mean(a[m] ** 2))


def main(png_dir, linear):
    wet3 = native(rdmds(HR_GRID + 'hFacC')) > 0
    mg = np.array([rms_v(MG_RUN, c, wet3) for c, _ in CTRLS])
    hr = np.array([rms_v(HR_RUN, c, wet3) for c, _ in CTRLS])
    for (c, _), a, b in zip(CTRLS, mg, hr):
        print('%-7s MG %.2e  HR %.2e  MG/HR %.2f' % (c, a, b, a / b if b > 0 else np.nan))

    # v is in units of the prior sigma; show it as a percentage of sigma
    mg, hr = 100 * mg, 100 * hr
    labels = [lab for _, lab in CTRLS]
    x = np.arange(len(labels))
    w = 0.38
    fig, ax = plt.subplots(figsize=(10.0, 5.2))
    for off, vals, col, lab in ((-w / 2 - 0.01, mg, C_MG, 'Multi-grid'),
                                (w / 2 + 0.01, hr, C_HR, 'HR-only')):
        if linear:
            ax.bar(x + off, vals, width=w, color=col, label=lab, zorder=3)
        else:
            ax.bar(x + off, vals - FLOOR, width=w, bottom=FLOOR, color=col, label=lab, zorder=3)
    if linear:
        ax.set_ylim(0, 7)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d%%'))
    else:
        ax.set_yscale('log')
        ax.set_ylim(FLOOR, 20)
        ax.set_yticks([1e-3, 1e-2, 1e-1, 1, 10])
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_yticklabels(['0.001%', '0.01%', '0.1%', '1%', '10%'])
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=13, rotation=35, ha='right', rotation_mode='anchor')
    ax.set_ylabel('RMS adjustment\n(% of prior uncertainty)', fontsize=14)
    ax.tick_params(axis='y', labelsize=12)
    ax.tick_params(axis='x', length=0)
    ax.grid(axis='y', which='major', color='0.88', lw=0.8, zorder=0)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=13, frameon=False, loc='upper right')

    name = NAME + ('_linear' if linear else '')
    fig.savefig(OUTDIR + name + '.pdf', bbox_inches='tight')
    print('wrote', OUTDIR + name + '.pdf')
    if png_dir:
        fig.savefig(os.path.join(png_dir, name + '.png'), dpi=200, bbox_inches='tight')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--png-dir', default=None)
    ap.add_argument('--linear', action='store_true', help='linear y-axis (writes *_linear)')
    a = ap.parse_args()
    plt.rcParams.update({'pdf.fonttype': 42})
    main(a.png_dir, a.linear)
