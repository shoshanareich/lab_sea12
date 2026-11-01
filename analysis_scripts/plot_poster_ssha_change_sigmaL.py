#!/usr/bin/env python3
"""
Poster drafts, 10-day sigmaL runs: what the optimization changed, HR-only vs MG.

1. ssha_change_swaths_HR_vs_MG_sigmaL_iter10      (along the SWOT swaths)
       SWOT - iter 0  |  MG iter 10 - iter 0  |  HR iter 10 - iter 0
   Each field has its own sample mean removed first, the same offset obsfit_cost.F
   removes, so the mean dynamic topography cancels.  "Asked for" vs "delivered".

2. ssha_change_map_HR_vs_MG_sigmaL_iter10         (gridded, day-10 snapshot)
       MG Eta(iter 10) - Eta(iter 0)  |  HR same
   Only the first and last Eta snapshots are written (dumpFreq = 0), so this is
   the end of the window, not a 10-day mean.  Wet-area mean removed.

3. ctrl_sst_adjustment_HR_vs_MG_sigmaL_iter10
       initial SST adjustment applied at 1/12 deg:  MG  |  HR
   xx_theta.effective is the full initial temperature, not the perturbation, so
   the adjustment is effective(iter 10) - effective(iter 0).

Map style matches plot_poster_grad_sst_sigmaL.py.

Usage:  python3 plot_poster_ssha_change_sigmaL.py [--png-dir DIR] [--aspect A]
"""

import argparse
import os
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from plot_poster_maps_sigmaL import (rdmds, native, load_mod, HR_RUN, MG_RUN, HR_GRID,
                                     OBS_FILE, OUTDIR, EXTENT)

ITER = 10
PC = ccrs.PlateCarree()
NE_LAND = cfeature.NaturalEarthFeature('physical', 'land', '50m')
VMAX_SSH = 0.4     # m
VMAX_SST = 0.6     # deg C
MSIZE = 1.5


def save(fig, name, png_dir):
    fig.savefig(OUTDIR + name + '.pdf', bbox_inches='tight')
    print('wrote', OUTDIR + name + '.pdf')
    if png_dir:
        fig.savefig(os.path.join(png_dir, name + '.png'), dpi=200, bbox_inches='tight')


def map_axes(ncol, aspect, width):
    fig, axs = plt.subplots(1, ncol, figsize=(width, 5.6), subplot_kw=dict(projection=PC),
                            gridspec_kw=dict(wspace=0.05))
    for ax in axs:
        ax.set_extent(EXTENT, crs=PC)
        ax.set_aspect(aspect)
    return fig, axs


def finish(fig, axs, titles, im, label):
    for ax, t in zip(axs, titles):
        ax.add_feature(NE_LAND, facecolor='lightgray', edgecolor='k', linewidth=0.8)
        gl = ax.gridlines(crs=PC, draw_labels=True, xlocs=range(-70, -40, 10),
                          ylocs=range(50, 80, 10), color='0.6', lw=0.6, alpha=0.7,
                          x_inline=False, y_inline=False, rotate_labels=False)
        gl.top_labels = gl.right_labels = False
        if ax is not axs[0]:
            gl.left_labels = False
        gl.xlabel_style = gl.ylabel_style = dict(size=11)
        ax.set_title(t, fontsize=16)
    cb = fig.colorbar(im, ax=axs, orientation='vertical', fraction=0.03, pad=0.02,
                      aspect=20, extend='both')
    cb.set_label(label, fontsize=14)
    cb.ax.tick_params(labelsize=12)


def swaths(png_dir, aspect):
    ds = xr.open_dataset(OBS_FILE)
    obs = ds.obs_val.values
    lon, lat, day = ds.sample_lon.values, ds.sample_lat.values, ds.obs_YYYYMMDD.values
    m0 = load_mod(MG_RUN, 0)
    mods = {'MG': load_mod(MG_RUN, ITER), 'HR': load_mod(HR_RUN, ITER)}
    use = (obs > -9999.0) & np.isfinite(m0)
    for m in mods.values():
        use &= np.isfinite(m)
    print('samples:', use.sum(), ' dates', day[use].min(), '-', day[use].max())

    anom = lambda a: a[use] - a[use].mean()
    fields = [anom(obs) - anom(m0)] + [anom(mods[k]) - anom(m0) for k in ('MG', 'HR')]
    for k, f in zip(('innovation', 'MG change', 'HR change'), fields):
        print('%-11s rms %.1f cm   corr with innovation %.3f'
              % (k, 100 * f.std(), np.corrcoef(f, fields[0])[0, 1]))

    fig, axs = map_axes(3, aspect, 15.5)
    for ax, f in zip(axs, fields):
        im = ax.scatter(lon[use], lat[use], c=f, s=MSIZE, cmap='RdBu_r', vmin=-VMAX_SSH,
                        vmax=VMAX_SSH, linewidths=0, transform=PC, rasterized=True)
    finish(fig, axs, ['SWOT − initial model', 'Multi-grid change', 'HR-only change'],
           im, 'SSHA (m)')
    save(fig, 'ssha_change_swaths_HR_vs_MG_sigmaL_iter%d' % ITER, png_dir)


def gridded(png_dir, aspect):
    hf = native(rdmds(HR_GRID + 'hFacC'))[0]
    rA = native(rdmds(HR_GRID + 'RAC'))
    xc, yc = native(rdmds(HR_GRID + 'XC')), native(rdmds(HR_GRID + 'YC'))
    wet = hf > 0
    nstep = 2881

    def deta(run):
        d = (native(rdmds(run + 'run_adhi_it%04d/Eta' % ITER, itrs=nstep))
             - native(rdmds(run + 'run_adhi_it0000/Eta', itrs=nstep)))
        d -= np.sum(d[wet] * rA[wet]) / np.sum(rA[wet])
        return np.where(wet, d, np.nan)

    fig, axs = map_axes(2, aspect, 11.0)
    for ax, run, lab in zip(axs, (MG_RUN, HR_RUN), ('MG', 'HR')):
        d = deta(run)
        print('dEta day 10 %s: rms %.1f cm' % (lab, 100 * np.nanstd(d)))
        im = ax.pcolormesh(xc, yc, d, cmap='RdBu_r', vmin=-VMAX_SSH, vmax=VMAX_SSH,
                           shading='nearest', transform=PC, rasterized=True)
    finish(fig, axs, ['Multi-grid change, day 10', 'HR-only change, day 10'],
           im, 'SSHA (m)')
    save(fig, 'ssha_change_map_HR_vs_MG_sigmaL_iter%d' % ITER, png_dir)


def sst_adjustment(png_dir, aspect):
    hf = native(rdmds(HR_GRID + 'hFacC'))[0]
    xc, yc = native(rdmds(HR_GRID + 'XC')), native(rdmds(HR_GRID + 'YC'))

    def dsst(run):
        get = lambda it: native(rdmds(run + 'run_adhi_it%04d/xx_theta.effective' % it,
                                      itrs=np.nan))[0]
        return np.where(hf > 0, get(ITER) - get(0), np.nan)

    fig, axs = map_axes(2, aspect, 11.0)
    for ax, run, lab in zip(axs, (MG_RUN, HR_RUN), ('MG', 'HR')):
        d = dsst(run)
        print('dSST %s: rms %.2f C, max |d| %.2f C'
              % (lab, np.nanstd(d), np.nanmax(np.abs(d))))
        im = ax.pcolormesh(xc, yc, d, cmap='RdBu_r', vmin=-VMAX_SST, vmax=VMAX_SST,
                           shading='nearest', transform=PC, rasterized=True)
    finish(fig, axs, ['Multi-grid', 'HR-only'], im, 'Initial SST adjustment (°C)')
    save(fig, 'ctrl_sst_adjustment_HR_vs_MG_sigmaL_iter%d' % ITER, png_dir)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--png-dir', default=None)
    ap.add_argument('--aspect', type=float, default=1.3)
    a = ap.parse_args()
    plt.rcParams.update({'pdf.fonttype': 42})
    swaths(a.png_dir, a.aspect)
    gridded(a.png_dir, a.aspect)
    sst_adjustment(a.png_dir, a.aspect)
