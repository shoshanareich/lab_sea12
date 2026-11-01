import glob
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal

root_dir = '/scratch/shoshi/labsea_MG_12/assim_swot_MG_10days_xxmult16_newCtrlW_v2/'
obs_file = '/home/shoshi/MITgcm_c69j/lab_sea12/input_binaries_hires/swot_obsfit_cycles_9thru11_labsea_L3v3.nc'
obs_name = 'swot_obsfit_cycles_9thru11_labsea_L3v3'
OUT = '/tmp/claude-62959/-home-shoshi-MITgcm-c69j/5a06a804-938c-4714-9be4-12f41d8364b3/scratchpad/'

ds = xr.open_dataset(obs_file); nobs = ds.dims['iOBS']
obs = ds.obs_val.values; line_id = ds.line_id.values

def load_equi(run_dir, nobs, obs_name=obs_name):
    mv = np.full(nobs, np.nan); mm = np.zeros(nobs)
    for fn in sorted(glob.glob(run_dir + obs_name + '.*.equi.nc')):
        d = xr.open_dataset(fn); i = d.sample_ind_glob.values.astype(int) - 1
        ok = (i >= 0) & (i < nobs)
        mv[i[ok]] = d.mod_val.values[ok]; mm[i[ok]] = d.mod_mask.values[ok]; d.close()
    return mv, mm

mv0, mm0 = load_equi(root_dir + 'run_adhi_it0000/', nobs)
mv5, mm5 = load_equi(root_dir + 'run_adhi_it0005/', nobs)
usable = (mm0 == 1) & (mm5 == 1) & np.isfinite(mv0) & np.isfinite(mv5) & np.isfinite(obs)

gr = pd.DataFrame({'cycle': ds.cycle.values, 'pass_no': ds.pass_no.values,
                   'track_id': ds.track_id.values}).groupby(['cycle','pass_no','track_id']).groups

MIN = 30
old_tracks, new_segs, gaps, gaps_per_track = [], [], [], []
for combo, ind in gr.items():
    i = np.asarray(ind)
    t0 = mv0[i]*mm0[i]; t5 = mv5[i]*mm5[i]
    if len(i) >= MIN and not np.isnan(obs[i]).any() and not np.isnan(t0).any() and not np.isnan(t5).any():
        old_tracks.append(i)
        d = np.diff(line_id[i]); g = d[d != 1]
        gaps.extend(g); gaps_per_track.append(len(g))
    j = i[usable[i]]
    if len(j) < 2: continue
    for s in np.split(j, np.where(np.diff(line_id[j]) != 1)[0] + 1):
        if len(s) >= MIN: new_segs.append(s)
gaps = np.array(gaps); gaps_per_track = np.array(gaps_per_track)

dx = 10.0; f_nyq = 1/(2*dx)
F = np.logspace(np.log10(0.0005), np.log10(f_nyq), 120)
def psd(x):
    N = len(x)
    fr, P = signal.periodogram(signal.detrend(x), fs=1/dx, window='hann', scaling='density')
    y = np.exp(np.interp(F, fr[1:], np.log(P[1:]))); y[F < 1/(N*dx)] = np.nan
    return y
def stack(S, v):
    M = np.array([psd(v[s]) for s in S])
    n = np.sum(~np.isnan(M), 0)
    return np.nanmean(M, 0), 1.96*np.nanstd(M, 0)/np.sqrt(np.maximum(n, 1))

SERIES = [(obs, 'black', 2.0, 'Observations'),
          (mv0, 'crimson', 1.5, 'Model Iter 0'),
          (mv5, 'royalblue', 1.5, 'Model Iter 5')]

# ================= FIGURE 1: old vs new spectra =================
fig, axes = plt.subplots(1, 2, figsize=(14, 6.2), sharey=True)
for ax, S, ttl in [(axes[0], old_tracks, 'OLD: whole tracks (gaps ignored)'),
                   (axes[1], new_segs,  'NEW: gap-free segments')]:
    for v, c, lw, lab in SERIES:
        m, ci = stack(S, v)
        ax.loglog(F, m, color=c, lw=lw, label=lab)
        ax.fill_between(F, np.maximum(m-ci, 1e-10), m+ci, color=c, alpha=0.15)
    sig = 0.02*0.208
    ax.axhline(2*sig**2*dx, color='0.45', ls=':', lw=1.5, label='SWOT noise floor (4 mm)')
    ax.axvspan(1/60., f_nyq, color='0.85', alpha=0.5, zorder=0)
    ax.set_xlabel('Wavenumber (cycles / km)')
    ax.set_title('%s\n%d segments, %d points' % (ttl, len(S), sum(len(s) for s in S)), fontsize=11)
    ax.grid(True, which='both', ls='--', alpha=0.4)
    ax.set_xlim(4e-4, 7e-2); ax.set_ylim(1e-6, 2e1)
    ax2 = ax.secondary_xaxis('top', functions=(lambda x: 1/x, lambda x: 1/x))
    ax2.set_xlabel('Wavelength (km)')
    ax2.set_xticks([1000, 500, 200, 100, 50, 30, 20])
    ax2.set_xticklabels(['1000','500','200','100','50','30','20'])
axes[0].set_ylabel('PSD (m$^2$ / cycles/km)')
axes[0].legend(loc='lower left', fontsize=9)
# ratio of the two observed spectra, as an inset annotation
mo_old, _ = stack(old_tracks, obs); mo_new, _ = stack(new_segs, obs)
for lam, xy in [(50, None), (25, None)]:
    i = np.argmin(np.abs(1/F - lam))
    axes[0].annotate('%.1fx high\nat %d km' % (mo_old[i]/mo_new[i], lam),
                     xy=(F[i], mo_old[i]), xytext=(F[i]*0.55, mo_old[i]*22),
                     fontsize=8, color='0.25', ha='center',
                     arrowprops=dict(arrowstyle='->', color='0.45', lw=0.9))
fig.suptitle('Effect of splitting tracks at along-track gaps  '
             '(MG 10-day, assim_swot_MG_10days_xxmult16_newCtrlW_v2)', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig(OUT+'cmp_spectra.png', dpi=150)

# ================= FIGURE 2: histograms =================
fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
L_old = np.array([len(s) for s in old_tracks]); L_new = np.array([len(s) for s in new_segs])
ax = axes[0]
b = np.arange(0, 275, 10)
ax.hist(L_old, bins=b, color='crimson', alpha=0.55, label='OLD whole tracks (n=%d)' % len(L_old))
ax.hist(L_new, bins=b, color='royalblue', alpha=0.55, label='NEW gap-free segments (n=%d)' % len(L_new))
ax.axvline(MIN, color='0.3', ls='--', lw=1.2)
ax.text(MIN+4, ax.get_ylim()[1]*0.93, 'min_len = %d' % MIN, fontsize=8, color='0.3')
ax.set_xlabel('points per segment  (x10 km)'); ax.set_ylabel('count')
ax.set_title('Segment length\nmedian %d -> %d points' % (np.median(L_old), np.median(L_new)), fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.3)

ax = axes[1]
b = np.arange(1.5, gaps.max()+2, 1)
ax.hist(gaps*dx, bins=b*dx, color='crimson', alpha=0.75)
ax.set_yscale('log')
ax.set_xlabel('gap size (km) spanned as if 10 km')
ax.set_ylabel('count')
ax.set_title('Gaps inside the OLD tracks\n%d gaps, max %d km' % (len(gaps), gaps.max()*dx), fontsize=10)
ax.grid(alpha=0.3)

ax = axes[2]
b = np.arange(-0.5, gaps_per_track.max()+1.5, 1)
ax.hist(gaps_per_track, bins=b, color='crimson', alpha=0.75)
ax.set_yscale('log')
ax.set_xlabel('gaps per OLD track'); ax.set_ylabel('count')
ax.set_title('%.0f%% of OLD tracks contain >= 1 gap' % (100*(gaps_per_track > 0).mean()), fontsize=10)
ax.grid(alpha=0.3)
fig.suptitle('Where the difference comes from', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT+'cmp_hist.png', dpi=150)

print('old %d tracks / new %d segs; gaps %d; tracks with gap %.0f%%; max gap %d km'
      % (len(old_tracks), len(new_segs), len(gaps), 100*(gaps_per_track>0).mean(), gaps.max()*dx))
for lam in (100, 50, 25):
    i = np.argmin(np.abs(1/F-lam)); print('  obs PSD ratio old/new at %d km: %.2f' % (lam, mo_old[i]/mo_new[i]))
