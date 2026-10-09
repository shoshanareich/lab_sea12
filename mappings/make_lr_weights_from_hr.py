#!/usr/bin/env python3
"""
Build the LR control-prior weight files by RESTRICTING the HR ones.

Why restrict instead of regenerating from ASTE at low resolution:

  The HR forward applies   u = Sigma_HR^(1/2) C^(1/2) P v   with P = piecewise-
  constant injection, while the LR adjoint differentiates a map built from
  Sigma_LR^(1/2).  For the LR gradient to be a gradient of the HR cost we need

      Sigma_LR = R Sigma_HR P                                             (*)

  With sigma = 1/sqrt(weight) and R the wet-area-weighted block mean, (*) holds
  exactly when

      sigma_LR(c) = sum_block( hFacC*rA * sigma_HR ) / sum_block( hFacC*rA )

  i.e. the area-weighted mean of SIGMA (not of sigma^2, and not the block max
  that wgenarr3d_lo uses -- both fail (*); see the adjoint test at the bottom,
  which this script runs every time).

  Regenerating at LR instead would (a) not satisfy (*), and (b) smooth with
  sigma=2 in LR cells, about 74 km rather than the 18 km used at HR -- the same
  nominal parameter, a 4x different physical prior.

Run after regenerating the HR weights.  Output is (Nr, Ny_lr, Nx_lr) big-endian
float64, the layout the LR model's READ_REC_3D_RL expects.
"""

import os
import sys

import numpy as np

sys.path.append('/home/shoshi/MITgcm_c68r/MITgcm/utils/python/MITgcmutils')
from MITgcmutils import rdmds

# ============================== configuration ==============================
HR_GRID = '/scratch/shoshi/labsea_MG_12/grid_hires_cleanbathy'
LR_GRID = '/scratch/shoshi/labsea_MG_12/grid_lores_cleanbathy_v2'

HR_DIR = '/home/shoshi/MITgcm_obsfit/lab_sea12/input_weights_hires/'
LR_DIR = '/home/shoshi/MITgcm_c69j/lab_sea12/input_weights_lores/'

FAC = 4                          # HR cells per LR cell, per direction (2 for 1/6 deg)

FIELDS = {                       # field: (HR input, LR output)
    'T': ('T_wt_fromASTE_HR_removeboundary_v2.bin', 'T_wt_fromASTE_restrictHR_v1.bin'),
    'S': ('S_wt_fromASTE_HR_removeboundary_v2.bin', 'S_wt_fromASTE_restrictHR_v1.bin'),
}

# the file currently named by input_adlo/data.ctrl, for comparison only
CURRENT_LR = {'T': 'T_wt_fromASTE_removeboundary_v4.bin',
              'S': 'S_wt_fromASTE_removeboundary_v4.bin'}
# ===========================================================================


def native(a, dtype='f8'):
    """MITgcm binaries are big-endian; keep everything native inside numpy."""
    return np.asarray(np.asarray(a), dtype=dtype)


hf_hr = native(rdmds(HR_GRID + '/hFacC'))
rA_hr = native(rdmds(HR_GRID + '/RAC'))
hf_lr = native(rdmds(LR_GRID + '/hFacC'))

NZ, NYH, NXH = hf_hr.shape
NY, NX = NYH // FAC, NXH // FAC
assert hf_lr.shape == (NZ, NY, NX), (hf_hr.shape, hf_lr.shape, FAC)
wet_lr = hf_lr > 0

# the LR bathymetry must be the coarsening of the HR one, or (*) is meaningless
coarse = (hf_hr > 0).reshape(NZ, NY, FAC, NX, FAC).any(axis=(2, 4))
assert (coarse == wet_lr).all(), 'LR hFacC is not the %dx%d coarsening of HR hFacC' % (FAC, FAC)
print('grids: HR %s -> LR %s, factor %d;  LR wet %d;  bathymetry consistent'
      % (hf_hr.shape, hf_lr.shape, FAC, int(wet_lr.sum())))


def restrict(w_hr):
    """Area-weighted block mean of sigma.  Returns (w_lr, sigma_lr, area, sigma_hr)."""
    use = (hf_hr > 0) & (w_hr > 0)                     # HR points that carry a control
    A = hf_hr * rA_hr[None, :, :] * use                # wet area, zero where unusable
    sigma_hr = np.where(use, 1.0 / np.sqrt(np.where(use, w_hr, 1.0)), 0.0)

    den = A.reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))
    num = (sigma_hr * A).reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))
    sigma_lr = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)

    w_lr = np.where(sigma_lr > 0, 1.0 / np.where(sigma_lr > 0, sigma_lr, 1.0) ** 2, 0.0)
    return w_lr, sigma_lr, A, sigma_hr


def adjoint_test(A, sigma_hr, sigma_lr, seed=0):
    """Check  R (Sigma_HR P v) == Sigma_LR v  on a random draw.  Should be ~1e-16."""
    rng = np.random.RandomState(seed)
    v = rng.randn(NZ, NY, NX)
    Pv = np.repeat(np.repeat(v, FAC, axis=-2), FAC, axis=-1)     # injection
    den = A.reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))
    lhs = np.where(den > 0,
                   (sigma_hr * Pv * A).reshape(NZ, NY, FAC, NX, FAC).sum(axis=(2, 4))
                   / np.where(den > 0, den, 1.0), 0.0)
    rhs = sigma_lr * v
    m = den > 0
    return np.abs(lhs[m] - rhs[m]).max() / np.abs(lhs[m]).max()


for field, (fin, fout) in FIELDS.items():
    src = HR_DIR + fin
    if not os.path.exists(src):
        print('\n=== %s: skipped, no %s' % (field, src))
        continue
    print('\n=== %s  <-  %s' % (field, fin))

    n = os.path.getsize(src)
    assert n == NZ * NYH * NXH * 8, 'expected %d bytes for >f8, got %d' % (NZ * NYH * NXH * 8, n)
    w_hr = np.fromfile(src, dtype='>f8').reshape(NZ, NYH, NXH)

    wet_hr = hf_hr > 0
    bad = int((wet_hr & (w_hr <= 0)).sum())
    print('  HR file: %d wet cells, %d with zero-or-negative weight%s'
          % (int(wet_hr.sum()), bad, '  <-- fix the HR file first' if bad else ''))

    w_lr, sigma_lr, A, sigma_hr = restrict(w_hr)

    err = adjoint_test(A, sigma_hr, sigma_lr)
    print('  adjoint test  R Sigma_HR P v == Sigma_LR v :  max rel err %.3e' % err)
    assert err < 1e-12, 'restriction is not the adjoint of injection'

    holes = int((wet_lr & (w_lr <= 0)).sum())
    print('  LR: %d wet cells, %d with no usable HR support (weight 0)'
          % (int(wet_lr.sum()), holes))
    print('  finite everywhere: %s   negatives: %d'
          % (bool(np.isfinite(w_lr).all()), int((w_lr < 0).sum())))

    s = sigma_lr[wet_lr & (sigma_lr > 0)]
    sh = sigma_hr[wet_hr & (sigma_hr > 0)]
    print('  sigma HR : median %.4g  p99 %.4g  max %.4g' % (np.median(sh),
          np.percentile(sh, 99), sh.max()))
    print('  sigma LR : median %.4g  p99 %.4g  max %.4g' % (np.median(s),
          np.percentile(s, 99), s.max()))

    cur = LR_DIR + CURRENT_LR[field]
    if os.path.exists(cur):
        wc = np.fromfile(cur, dtype='>f8').reshape(NZ, NY, NX)
        m = wet_lr & (wc > 0) & (w_lr > 0)
        r = (1 / np.sqrt(wc[m])) / (1 / np.sqrt(w_lr[m]))
        print('  vs %s on %d shared cells:' % (CURRENT_LR[field], int(m.sum())))
        print('     sigma_current / sigma_new : median %.4f  p95 %.3f  p99 %.3f  max %.1f'
              % (np.median(r), np.percentile(r, 95), np.percentile(r, 99), r.max()))
        print('     cells the current file has but the new one drops: %d'
              % int((wet_lr & (wc > 0) & (w_lr <= 0)).sum()))

    w_lr.astype('>f8').tofile(LR_DIR + fout)
    print('  wrote %s%s  (%d bytes)' % (LR_DIR, fout, os.path.getsize(LR_DIR + fout)))

print("""
Next: point input_adlo/data.ctrl at the new files

  xx_genarr3d_weight(1) = '%s',
  xx_genarr3d_weight(2) = '%s',

and start a FRESH run directory from iteration 0 with coldstart=.TRUE.  The
penalty is sum(v^2) with no weight, so changing sigma silently reinterprets any
existing ecco_ctrl_*.opt* and OPWARM.  The 12/9/11/16 mult_gen* factors in the
xxmult16 driver were measured under the old weights and want re-deriving too.
""" % (FIELDS['T'][1], FIELDS['S'][1]))
