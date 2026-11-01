#!/usr/bin/env python3
"""
Poster figure: initial-SST gradient, HR-only vs multi-grid adjoint, 10-day
sigmaL runs, iteration 0.

    gradient_consistency_SST_HR_vs_MG_sigmaL_iter0

adxx_theta level 0 for both runs, unscaled, on one shared colour scale.  The
HR-only 1/12 deg gradient is restricted to 1/3 deg with the wet-area-weighted
block mean, as in plot_poster_maps_sigmaL.py.  Plate carree: the model
domain is a lon/lat box, so it fills the frame exactly.  Land is Natural Earth 50m only,
so model-closed basins (Hudson Strait, Gulf of St. Lawrence) show as ocean.

Usage:  python3 plot_poster_grad_sst_sigmaL.py [--png-dir DIR] [--vmax V] [--aspect A]
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from plot_poster_maps_sigmaL import (rdmds, native, restrict2d, HR_RUN, MG_RUN,
                                     HR_GRID, LR_GRID, OUTDIR, NX, NY, FAC, LANDC,
                                     EXTENT)

NAME = 'gradient_consistency_SST_HR_vs_MG_sigmaL_iter0'
PC = ccrs.PlateCarree()
PROJ = PC
OCEANC = 'white'     # ocean outside the model's wet cells


def main(png_dir, vmax, aspect):
    hf_hr = native(rdmds(HR_GRID + 'hFacC'))[0]
    rA_hr = native(rdmds(HR_GRID + 'RAC'))
    hf_lr = native(rdmds(LR_GRID + 'hFacC'))[0]
    xl, yl = native(rdmds(LR_GRID + 'XC')), native(rdmds(LR_GRID + 'YC'))
    W_hr = hf_hr * rA_hr
    W_lr = W_hr.reshape(NY, FAC, NX, FAC).sum(axis=(1, 3))
    wet = (hf_lr > 0) & (W_lr > 0)

    g_hr = native(rdmds(HR_RUN + 'run_adhi_it0000/adxx_theta', itrs=np.nan))[0]
    g_mg = native(rdmds(MG_RUN + 'run_adlo_it0000/adxx_theta', itrs=np.nan))[0]
    g_hr = np.where(wet, restrict2d(g_hr, W_hr, W_lr), np.nan)
    g_mg = np.where(wet, g_mg, np.nan)
    for k, g in (('HR', g_hr), ('MG', g_mg)):
        print('%s  |g| 50/99/100th pct: %s' % (k, np.nanpercentile(np.abs(g), [50, 99, 100])))

    ne_land = cfeature.NaturalEarthFeature('physical', 'land', '50m')
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 5.6), subplot_kw=dict(projection=PROJ),
                            gridspec_kw=dict(wspace=0.05))
    for ax, g, lab in ((axs[0], g_hr, 'HR-only adjoint'),
                       (axs[1], g_mg, 'Multi-grid adjoint')):
        ax.set_extent(EXTENT, crs=PC)
        ax.set_aspect(aspect)   # >1 stretches latitude: taller panels for the poster
        ax.set_facecolor(OCEANC)
        im = ax.pcolormesh(xl, yl, g, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                           shading='nearest', transform=PC, rasterized=True)
        ax.add_feature(ne_land, facecolor='lightgray', edgecolor='k', linewidth=0.8)
        gl = ax.gridlines(crs=PC, draw_labels=True, xlocs=range(-70, -40, 10),
                          ylocs=range(50, 80, 10), color='0.6', lw=0.6, alpha=0.7,
                          x_inline=False, y_inline=False, rotate_labels=False)
        gl.top_labels = gl.right_labels = False
        if ax is axs[1]:
            gl.left_labels = False
        gl.xlabel_style = gl.ylabel_style = dict(size=11)
        ax.set_title(lab, fontsize=16)
    cb = fig.colorbar(im, ax=axs, orientation='vertical', fraction=0.03, pad=0.02,
                      aspect=20, extend='both')
    cb.set_label(r'$\partial J\,/\,\partial\,\mathrm{SST}$', fontsize=14)
    cb.ax.tick_params(labelsize=12)

    fig.savefig(OUTDIR + NAME + '.pdf', bbox_inches='tight')
    print('wrote', OUTDIR + NAME + '.pdf')
    if png_dir:
        fig.savefig(os.path.join(png_dir, NAME + '.png'), dpi=200, bbox_inches='tight')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--png-dir', default=None)
    ap.add_argument('--vmax', type=float, default=30.0)
    ap.add_argument('--aspect', type=float, default=1.3,
                    help='height of 1 deg lat / 1 deg lon in the panels')
    a = ap.parse_args()
    plt.rcParams.update({'pdf.fonttype': 42})
    main(a.png_dir, a.vmax, a.aspect)
