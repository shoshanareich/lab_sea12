# Multi-grid adjoint SWOT assimilation in the Labrador Sea

Assimilation of **SWOT** sea-surface-height observations (plus Argo T/S profiles) into a
regional **Labrador Sea** configuration of **MITgcm** (checkpoint **c69j**), using the ECCO
framework (`pkg/ctrl`, `pkg/ecco`, `pkg/obsfit`, `pkg/profiles`) and a TAF-generated adjoint.

SWOT resolves very small scales, so the forward model has to run at high resolution. An
adjoint at that resolution is expensive, so this repo uses a **multi-grid (MG) adjoint**:

| | role | resolution | grid | `deltaT` | procs |
|---|---|---|---|---|---|
| **HR** (`adhi`) | forward model | 1/12° | 480 × 384 × 21 | 300 s | 103 |
| **LR** (`adlo`) | adjoint model | 1/3° | 120 × 96 × 21 | 1200 s | 16 |

The two are joined by MITgcm's **divided adjoint (DivA)** plus a local modification,
`MULTISCALE_COUPLING_TAPES` that coarsens the HR checkpoints (tapes) to the LR. The HR forward run writes its checkpoint tapes already
**coarsened** by `MSfacX × MSfacY = 4 × 4` blocks, and the LR adjoint reads them as its
reference trajectory. The optimizer and the control vector live at LR. Control adjustments
are prolongated back to HR for the next outer iteration.

An experimental 1/6° adjoint level (`adlo6`, `MSfac = 2`) is also kept here.

---

## One outer iteration

```
HR forward (1/12°, mitgcmuv_ad, DivA forward leg)
    → coarsened tapes/, cost, model equivalents at observation points
        ↓ mappings/make_cost_cp_v2.py          block-average ×4 → LR  (gencost misfits, xx/adxx)
        ↓ mappings/make_obsfit_lr_tiles.py     re-tile SWOT model equivalents → LR tiling
        ↓ mappings/make_prof_lr_tiles.py       re-tile Argo model equivalents → LR tiling
LR adjoint (1/3°, mitgcmuv_ad, DivA reverse sweep, invoked twice)
    → gradient at LR
OPTIM (optim.x / m1qn3, entirely at LR)
    → ecco_ctrl_MIT_CE_000.opt<iter+1>
pack/unpack (LR build, doMainUnpack=.TRUE.)
    → xx_* at LR
        ↓ mappings/interp_xx_lores_to_hires_itX_v2.py   prolongation LR → HR
xx_hires/  →  read by the next iteration's HR forward
```

The gradient is never computed at HR, so `adxx_*` in an HR run directory is identically
zero.

**The prolongation and restriction must be mutually consistent.** For the LR gradient to be a
gradient of the HR cost, the LR→HR prolongation has to match the LR control
parameterization (piecewise-constant **injection**), and the HR→LR restriction has to be its
adjoint (a **wet-area-weighted block mean**, `Σ hFacC·rA·w / Σ hFacC·rA`). Both are the
defaults:

| script | argument | values |
|---|---|---|
| `interp_xx_lores_to_hires_itX_v2.py` | argv 7 | `injection` (default) \| `linear` |
| | argv 8 | HR grid dir (default `grid_hires_cleanbathy`) |
| `make_cost_cp_v2.py` | argv 9 | `area` (default) \| `mean` |
| | argv 10 | HR grid dir (default `grid_hires_cleanbathy`) |

The old `linear`/`mean` pair is not adjoint-consistent and breaks the m1qn3 line search. It is
kept only for reproducing pre-September-2026 runs.

---

## Repository layout

This repo sits next to a **pristine** upstream MITgcm checkout:

```
MITgcm_c69j/
├── MITgcm/        upstream MITgcm c69j — never edited
└── lab_sea12/     this repo
```

All customization goes through `-mods` override directories here. If an upstream file needs
changing, copy it into the relevant `code_*` directory.

| path | contents |
|---|---|
| `code_adhi_2lev_seaice_update/` | HR (1/12°) adjoint-build overrides, used as the DivA forward leg |
| `code_adlo_2lev_seaice_update/` | LR (1/3°) adjoint-build overrides, which run the reverse sweep |
| `code_adlo6_2lev_seaice_update/` | experimental 1/6° adjoint level |
| `input_adhi/`, `input_adlo/`, `input_adlo6/` | namelists (`data*`, `eedata`) for each level |
| `mappings/` | offline HR↔LR Python transfers (see below) |
| `analysis_scripts/` | plotting and diagnostics: cost curves, gradients, control increments, SWOT spectra |
| `run_*.sh` | optimization drivers, one per experiment |
| `link_hires*.sh`, `link_lores*.sh`, `link_packunpack*.sh` | stage each run directory |
| `submit*.sh` | scheduler wrappers and standalone job files |
| `machine_config.sh`, `mappings/machine_config.py` | **the only** places machine-specific paths live |
| `check_staging.sh` | dry-run checker for the link scripts |

`.gitignore` ignores everything at the top level by default and whitelists only the items
above. Build directories, run output, binaries, forcing, weights and grids are deliberately
untracked.

### Where the multi-grid machinery lives

All of it is guarded by `#ifdef MULTISCALE_COUPLING_TAPES`. None of it exists upstream.

- `mdsio_write_tape.F` (HR) writes tapes at `Nx*Ny/MSfacX/MSfacY`, block-averaging and
  skipping fill values.
- `mdsio_read_tape.F` (LR) reads the coarsened tapes.
- `AUTODIFF_PARAMS.h` and `autodiff_readparms.F` handle `MSfacX`/`MSfacY` from
  `data.autodiff`.
- `the_main_loop.F` adds `EXCH_*` calls before the inner checkpoint loop, because coarsened
  tapes arrive without valid halos (marked `CMM`).
- `ctrl_*`, `cost_generic.F`, `ecco_toolbox.F`, `packages_init_variables.F` and
  `convective_adjustment_ini.F` make the control and cost handling multiscale-aware.

**The shared files must stay in sync between `code_adhi_*` and `code_adlo_*`.** A tape-format
change on one side silently breaks the other, so diff them after editing.

### `mappings/`

| script | purpose |
|---|---|
| `make_cost_cp_v2.py` | restricts gridded `pkg/ecco` misfits and `xx_*`/`adxx_*`/`m_*`/`adm_*` from HR to LR |
| `make_obsfit_lr_tiles.py` | re-assigns SWOT samples (`pkg/obsfit`) from HR to LR tiles; moves model equivalents, not observations |
| `make_prof_lr_tiles.py` | same for Argo profiles (`pkg/profiles`) |
| `interp_xx_lores_to_hires_itX_v2.py` | prolongs the new LR control vector to HR |
| `make_lr_weights_from_hr.py` | builds LR control prior weights by restricting the HR ones; self-checks with a dot-product test |
| `map_hr_lr_tiles.py`, `file_utils.py`, `read_write.py` | helpers |
| `*_itX.py`, `make_cost_cp.py` (no `_v2`) | older versions, kept for reference |

The scripts run in conda env `py38` and hardcode `nx=120, ny=96, nz=21, factor=4` at the
top. A resolution change means editing those values.

---

## Configuration summary

| | HR | LR |
|---|---|---|
| `nTimeSteps` (63 days) | 18144 | 4536 |
| `nchklev_1 × nchklev_2` (`tamc.h`) | 48 × 378 | 12 × 378 |
| `MSfacX`, `MSfacY` | 4, 4 | 1, 1 |
| tiling (`SIZE.h` + `data.exch2`) | 32×32 tiles, 103 procs | 5×4 tiles − 4 blanks, 16 procs |
| packages | `gfd kpp exf seaice salt_plume diagnostics exch2 profiles obsfit ecco adjoint smooth` | same but `oceanic` instead of `gfd kpp`, **no `smooth`** |

- `nchklev_2` must match between HR and LR so the DivA chunk boundaries line up.
- `SIZE.h` and `data.exch2` must describe the same tiling. Alternates (`SIZE.h_48x48` /
  `data.exch2_48x48`, …) are kept side by side, so copy **both** when switching.
- `pkg/smooth` (the `WC01` correlation operator) runs **only at HR**, on purpose. Do not add it
  to the LR build for MG runs. It **is** needed in a single-resolution low-res optimization,
  otherwise the 3D initial-condition controls are unregularized.

### Cost and controls

- **SWOT** through `pkg/obsfit` (`data.obsfit`), **Argo** through `pkg/profiles`
  (`data.profiles`). The `pkg/ecco` gencost SSH terms are present but at weight 0.
- **Controls** (`data.ctrl`): daily `gentim2d` atmospheric adjustments (`xx_atemp`,
  `xx_precip`, `xx_swdown`, `xx_lwdown`, `xx_uwind`, `xx_vwind`, `xx_aqh`), plus `genarr3d`
  initial-condition controls (`xx_theta`, `xx_salt`, `xx_diffkr`) where enabled.
- **Prior weights** (`weight = 1/σ²`):
  - 3D T/S come from ASTE, interpolated to HR.
  - Atmospheric weights are the std of deseasonalized 2024 JRA-3Q daily means.
  - In both cases the LR weights are the **restriction of the HR weights**, not an
    independent regeneration.
- Suffixed namelists (`data_atmONLY`, `data.ctrl_atmONLY`, …) are swapped in by the matching
  `link_*` script. Check which link script a driver calls before reading a namelist.
- The `mult_*` factors in `data.obsfit` / `data.profiles` / `data.ecco` decide which
  observations actually enter the cost. **`input_adhi/data.obsfit` ships with SWOT off**
  (`mult_obsfit(1) = 0.0`), so drivers must set it.

---

## Building

From a build directory alongside the code directory:

```bash
cd lab_sea12/build_adhi_2lev_seaice_update_mpi
../../MITgcm/tools/genmake2 -mods ../code_adhi_2lev_seaice_update/ \
    -optfile=<optfile from machine_config.sh> -mpi -rd ../../MITgcm/ \
    -adof ../../MITgcm/tools/adjoint_options/adjoint_default
make depend
make adall          # → mitgcmuv_ad
```

- Adjoint builds use DivA (`USE_DIVA=1`) and run TAF on the remote FastOpt server, so they
  need network access. If one fails, check `taf_ad.log`.
- Forward-only builds (`code_hi_fwd`, `code_lo_fwd`, used for spinups and grids) use the same
  command without `-adof`, then `make all`.
- Rebuild after changing any of `SIZE.h`, `tamc.h`, `*_OPTIONS.h`, `CTRL_SIZE.h`,
  `OBSFIT_SIZE.h`, `PROFILES_SIZE.h` or `DIAGNOSTICS_SIZE.h`.
- After a long gap, run `ldd mitgcmuv_ad`. Cluster netcdf module updates break old binaries;
  rebuild rather than patching `LD_LIBRARY_PATH`.

---

## Running

### Machines

The drivers are machine-independent. Paths, modules, scheduler, MPI launcher and conda
activation all come from `machine_config.sh` (shell) and `mappings/machine_config.py`
(Python). Keep the two in sync. The machine is auto-detected (sverdrup, TACC Stampede3, or
NASA Pleiades/pfe), or forced with `export LABSEA_MACHINE=sverdrup|stampede3|pfe`.

### Submitting

```bash
./submit.sh run_optimization.sh        # detects machine → sbatch or qsub
```

This dispatches to `submit_slurm_sverdrup.sh`, `submit_slurm_stampede3.sh` or
`submit_pbs_pfe.sh`. You can also call a wrapper directly:

```bash
# sverdrup / stampede3
sbatch --export=ALL,DRIVER=run_optimization.sh submit_slurm_sverdrup.sh
# pfe
qsub -v DRIVER=run_optimization.sh submit_pbs_pfe.sh
```

Other useful forms:

```bash
./submit.sh                                         # no driver → list available ones
./submit.sh -n run_optimization.sh                  # dry run: print the command only
./submit.sh run_swot_optimization.sh -t 12:00:00    # extra args go to sbatch/qsub
LABSEA_MACHINE=pfe ./submit.sh -n run_optimization.sh   # test another machine's config
./run_optimization.sh                               # directly, inside an existing allocation
./submit.sh run_selftest.sh -t 00:05:00             # ~1 min end-to-end check of the submit chain
```

`run_optimization.sh` (Argo) is the cleanest reference for the loop.
`run_swot_optimization.sh` is the SWOT version. The other `run_*` drivers (`run_LR_opt.sh`,
`run_LR6*.sh`, `run_ecco_optimization.sh`, `run_summer*_optimization.sh`) differ mainly in
`scratchdir`, link-script variant and `data.optim`.

Some `submit_*.sh` files are **standalone job files** with hardcoded paths. A SLURM job file
cannot source `machine_config.sh`, because it runs from the spool directory. Those files are
not portable, and most are untracked.

### Checking staging without running

```bash
./check_staging.sh                    # check every link script against this machine
./check_staging.sh -s link_lores.sh   # one script
./check_staging.sh -v                 # list every source, not just problems
./check_staging.sh -m pfe             # path formation only for another machine
```

This runs the real `link_*.sh` control flow in a sandbox with file operations stubbed out, and
reports any missing source. It needs no allocation.

### The divided-adjoint restart

The LR adjoint is launched twice. The first run does one DivA level and writes
`divided.ctrl` as `377 376`. The driver then sets the end level to 0
(`sed -i 's/376/0/g' divided.ctrl`), and the second run finishes the sweep. **The `376`
literal is tied to `nchklev_2 = 378`.** If `nchklev_2` changes, update it in every driver,
because a stale literal fails silently.

---

## Gotchas worth knowing up front

- **`MSfacX = MSfacY = 1` for any single-resolution adjoint.** `input_adhi/data.autodiff`
  ships `4,4`. With coarsened tapes, a same-resolution adjoint dies in `mds_read_tape` with
  `attempt to access non-existent record`.
- **The operator defaults flipped in September 2026** (`linear`/`mean` →
  `injection`/`area`). Drivers that pass no operator arguments will not reproduce their
  pre-September output.
- **m1qn3 warm start:** use `coldstart=.TRUE.` only on the first iteration, and don't delete
  `OPWARM`. A non-monotone cost history is usually the line search probing.
- **Changing prior weights reinterprets the control vector.** Start a fresh run directory at
  iteration 0.
- Binaries are a mix of 32- and 64-bit (`readBinaryPrec=32`, but `*_64b.bin` files exist).
  Check `.meta` before reading raw files.
- `useSingleCpuIO=.TRUE.`, so output files are global, not per-tile.
