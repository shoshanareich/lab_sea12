#!/usr/bin/env python3
"""
Poster spatial figures, 10-day sigmaL runs.

1. ssh_misfit_maps_HR_vs_MG_sigmaL_iter10
       row 1  SSH anomaly:      SWOT | multi-grid iter 10 | HR-only iter 10
       row 2  normalized misfit z:  iter 0 | multi-grid iter 10 | HR-only iter 10
   z = (mod - obs - offset) / obs_uncert, offset = <mod> - <obs> as in
   obsfit_cost.F (see plot_optim_results_HR_vs_MG.py).  obs_uncert is the
   sigmaL value, so these z are not comparable to the ctrlWyearly figures.
   The 3 samples the sigmaL file sets to -9999 are dropped.

2. gradient_consistency_HR_vs_MG_sigmaL_iter0
       rows   adxx_theta level 0 (initial SST)  |  adxx_atemp summed over days
       cols   HR-only 1/12 deg adjoint, restricted to 1/3 deg  |  MG 1/3 deg adjoint
   Iteration 0 only: both runs share the identical state there, so the panels
   differ only in the adjoint used.  Restriction is the wet-area-weighted block
   mean (plot_grad_hr_vs_mg.py).  Each panel is scaled by its own 99th
   percentile of |g| -- the two gradients are with respect to controls on
   different grids, so the comparison is of pattern; the pattern correlation
   over wet 1/3 deg cells is printed and drawn on the figure.

Usage:  python3 plot_poster_maps_sigmaL.py [--png-dir DIR]
"""

import argparse
import os
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds

RUNSROOT = '/scratch/shoshi/labsea_MG_12/'
HR_RUN   = RUNSROOT + 'assim_swot_HR_10days_sigmaL/'
MG_RUN   = RUNSROOT + 'assim_swot_MG_10days_sigmaL/'
HR_GRID  = RUNSROOT + 'grid_hires_cleanbathy/'
LR_GRID  = RUNSROOT + 'grid_lores_cleanbathy_v2/'
OBSNAME  = 'swot_obsfit_cycles_9thru11_labsea_L3v3_sigmaL'
OBS_FILE = '/home/shoshi/MITgcm_c69j/lab_sea12/input_binaries_hires/' + OBSNAME + '.nc'
OUTDIR   = '/home/shoshi/MITgcm_c69j/lab_sea12/figs/'

ITER  = 10
NX, NY, FAC = 120, 96, 4
VMAX_SSH, VMAX_Z = 0.6, 5.0
MSIZE = 1.2
LANDC = '#bdbdbd'
EXTENT = (280, 320, 46, 78)


def native(a):
    return np.asarray(np.asarray(a), dtype='f8')


def save(fig, name, png_dir):
    fig.savefig(OUTDIR + name + '.pdf', bbox_inches='tight')
    print('wrote', OUTDIR + name + '.pdf')
    if png_dir:
        fig.savefig(os.path.join(png_dir, name + '.png'), dpi=200, bbox_inches='tight')


def land_hr(ax, hf, xc, yc):
    ax.pcolormesh(xc, yc, np.where(hf > 0, np.nan, 1.0), cmap=matplotlib.colors.ListedColormap([LANDC]),
                  shading='nearest', rasterized=True)
    ax.set_xlim(*EXTENT[:2]); ax.set_ylim(*EXTENT[2:])
    ax.set_aspect(1 / np.cos(np.deg2rad(62)))


# ------------------------------------------------------------- figure 1
def load_mod(run, it):
    fn = run + 'run_adhi_it%04d/' % it + OBSNAME + '.equi.nc'
    with xr.open_dataset(fn) as d:
        return d.mod_val.values.copy()


def ssh_misfit(png_dir):
    hf = native(rdmds(HR_GRID + 'hFacC'))[0]
    xc, yc = native(rdmds(HR_GRID + 'XC')), native(rdmds(HR_GRID + 'YC'))

    ds = xr.open_dataset(OBS_FILE)
    obs, unc = ds.obs_val.values, ds.obs_uncert.values
    lon, lat = ds.sample_lon.values, ds.sample_lat.values

    mods = {'it0': load_mod(MG_RUN, 0), 'MG': load_mod(MG_RUN, ITER),
            'HR': load_mod(HR_RUN, ITER)}
    use = (obs > -9999.0) & (unc > 0)
    for m in mods.values():
        use &= np.isfinite(m)
    print('samples:', use.sum())

    o = obs[use]; u = unc[use]
    ssh = {'SWOT': o - o.mean()}
    z = {}
    for k, m in mods.items():
        m = m[use]
        ssh[k] = m - m.mean()
        z[k] = (m - o - (m.mean() - o.mean())) / u
        print('%-4s <z^2> = %.2f   rms misfit = %.1f cm'
              % (k, np.mean(z[k] ** 2), 100 * np.sqrt(np.mean((z[k] * u) ** 2))))
    sx, sy = lon[use], lat[use]

    fig, axs = plt.subplots(2, 3, figsize=(11.5, 11.2), sharex=True, sharey=True,
                            gridspec_kw=dict(wspace=0.04, hspace=0.12))
    top = [('SWOT', 'SWOT SSH'), ('MG', 'Multi-grid, iter 10'), ('HR', 'HR-only, iter 10')]
    bot = [('it0', 'Initial model, iter 0'), ('MG', 'Multi-grid, iter 10'),
           ('HR', 'HR-only, iter 10')]
    for ax, (k, t) in zip(axs[0], top):
        land_hr(ax, hf, xc, yc)
        im1 = ax.scatter(sx, sy, c=ssh[k], s=MSIZE, cmap='RdBu_r', vmin=-VMAX_SSH,
                         vmax=VMAX_SSH, linewidths=0, rasterized=True)
        ax.set_title(t, fontsize=17)
    for ax, (k, t) in zip(axs[1], bot):
        land_hr(ax, hf, xc, yc)
        im2 = ax.scatter(sx, sy, c=z[k], s=MSIZE, cmap='RdBu_r', vmin=-VMAX_Z,
                         vmax=VMAX_Z, linewidths=0, rasterized=True)
        ax.set_title(t, fontsize=17)
        ax.text(0.97, 0.97, r'$\langle z^2\rangle$ = %.1f' % np.mean(z[k] ** 2),
                transform=ax.transAxes, ha='right', va='top', fontsize=15,
                bbox=dict(fc='white', ec='none', alpha=0.85, pad=3))
    for ax in axs[:, 0]:
        ax.set_ylabel('Latitude (°N)', fontsize=14)
    for ax in axs[1]:
        ax.set_xlabel('Longitude (°E)', fontsize=14)
    for ax in axs.flat:
        ax.tick_params(labelsize=12)

    c1 = fig.colorbar(im1, ax=axs[0], fraction=0.025, pad=0.015, extend='both')
    c1.set_label('SSH anomaly (m)', fontsize=14)
    c2 = fig.colorbar(im2, ax=axs[1], fraction=0.025, pad=0.015, extend='both',
                      ticks=[-5, -3, -1, 0, 1, 3, 5])
    c2.set_label(r'Normalized misfit  $z$', fontsize=14)
    for c in (c1, c2):
        c.ax.tick_params(labelsize=12)
    save(fig, 'ssh_misfit_maps_HR_vs_MG_sigmaL_iter10', png_dir)


# ------------------------------------------------------------- figure 2
def restrict2d(g_hr, W_hr, W_lr):
    num = (np.nan_to_num(g_hr) * W_hr).reshape(NY, FAC, NX, FAC).sum(axis=(1, 3))
    return np.where(W_lr > 0, num / np.where(W_lr > 0, W_lr, 1.0), np.nan)


def gradients(png_dir):
    hf_hr = native(rdmds(HR_GRID + 'hFacC'))[0]
    rA_hr = native(rdmds(HR_GRID + 'RAC'))
    hf_lr = native(rdmds(LR_GRID + 'hFacC'))[0]
    xl, yl = native(rdmds(LR_GRID + 'XC')), native(rdmds(LR_GRID + 'YC'))
    W_hr = hf_hr * rA_hr
    W_lr = W_hr.reshape(NY, FAC, NX, FAC).sum(axis=(1, 3))
    wet = (hf_lr > 0) & (W_lr > 0)

    def get(run, sub, field):
        a = native(rdmds(run + sub + '/' + field, itrs=np.nan))
        if field == 'adxx_theta':
            return a[0]
        return a.sum(axis=0) if a.ndim == 3 else a

    rows = [('adxx_theta', 'Initial SST'),
            ('adxx_atemp', 'Air temperature\n(summed over 10 days)')]
    fig, axs = plt.subplots(2, 2, figsize=(8.6, 10.4), sharex=True, sharey=True,
                            gridspec_kw=dict(wspace=0.04, hspace=0.12))
    for r, (field, name) in enumerate(rows):
        g_hr = restrict2d(get(HR_RUN, 'run_adhi_it0000', field), W_hr, W_lr)
        g_mg = get(MG_RUN, 'run_adlo_it0000', field)
        g_mg = np.where(wet, g_mg, np.nan)
        g_hr = np.where(wet, g_hr, np.nan)
        ok = wet & np.isfinite(g_hr) & np.isfinite(g_mg)
        corr = np.corrcoef(g_hr[ok], g_mg[ok])[0, 1]
        print('%-11s pattern correlation HR(restricted) vs MG, iter 0: %.3f' % (field, corr))
        for c, (g, lab) in enumerate([(g_hr, 'HR-only adjoint\n(1/12°, restricted to 1/3°)'),
                                      (g_mg, 'Multi-grid adjoint\n(1/3°)')]):
            ax = axs[r, c]
            s = np.nanpercentile(np.abs(g), 99)
            ax.set_facecolor(LANDC)
            im = ax.pcolormesh(xl, yl, g / s, cmap='RdBu_r', vmin=-1, vmax=1,
                               shading='nearest', rasterized=True)
            ax.set_xlim(*EXTENT[:2]); ax.set_ylim(*EXTENT[2:])
            ax.set_aspect(1 / np.cos(np.deg2rad(62)))
            ax.set_title(lab if r == 0 else '', fontsize=15)
            ax.tick_params(labelsize=12)
        axs[r, 0].set_ylabel(name + '\n\nLatitude (°N)', fontsize=14)
        axs[r, 1].text(0.97, 0.97, 'pattern\ncorr. = %.2f' % corr,
                       transform=axs[r, 1].transAxes, ha='right', va='top',
                       fontsize=15, bbox=dict(fc='white', ec='none', alpha=0.85, pad=3))
    for ax in axs[1]:
        ax.set_xlabel('Longitude (°E)', fontsize=14)
    cb = fig.colorbar(im, ax=axs, fraction=0.03, pad=0.02, extend='both',
                      ticks=[-1, -0.5, 0, 0.5, 1])
    cb.set_label(r'$\partial J/\partial v$, scaled by 99th percentile of $|g|$', fontsize=14)
    cb.ax.tick_params(labelsize=12)
    save(fig, 'gradient_consistency_HR_vs_MG_sigmaL_iter0', png_dir)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--png-dir', default=None)
    a = ap.parse_args()
    plt.rcParams.update({'pdf.fonttype': 42})
    ssh_misfit(a.png_dir)
    gradients(a.png_dir)
