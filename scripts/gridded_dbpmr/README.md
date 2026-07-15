# Gridded DBPM (dbpmr engine) — spatiotemporal, vertically biomass-weighted

A gridded extension of the LME/FAO DBPM calibration that runs **every 1° ocean cell** as an
independent 0-D `dbpmr` column (compiled-C size-spectrum engine), coupled only through an annual
**spatial-gravity** re-allocation of the region's fishing effort. It is the **engine-swap** analog of
the Python gridded pipeline (`07_setup_gridded_DBPM.py` / `08_run_dbpm_gridded.py`), using the `dbpmr`
solver and **per-cell, per-timestep, vertically BIOMASS-WEIGHTED** forcing.

Full spec + assumptions: [`GRIDDED_WORKFLOW.md`](GRIDDED_WORKFLOW.md).

## Key methodology (what's new vs the 0-D calibration)

1. **Vertically biomass-weighted forcing (no climatology, no depth-weighting).** Predators feed through
   the water column where the food is, so every plankton/temperature input is a phytoplankton-carbon
   weighted vertical mean over 0–200 m, per cell **per timestep**:
   `⟨X⟩ᵢ(t) = Σ_z Xᵢ(z,t)·phycᵢ(z,t)·Δz / Σ_z phycᵢ(z,t)·Δz`
   → intercept(t), slope(t) via `GetPPIntSlope`, and experienced temperature `T_exp=⟨thetao⟩`.
   Implemented as `integrating_phyto(..., weighting="biomass")` in `../useful_functions.py`.
2. **LME-centering** so the biomass-weighted regional aggregate reproduces the exact 0-D calibration
   input series (`intercept+dint`, `slope+dslope`, `tos+texp_offset`) — the calibrated q is preserved.
3. **Gridded q refit.** Resolving spatial+seasonal productivity heterogeneity inflates fish biomass
   (Jensen's inequality on a convex trophic response — FAO 58 ~26× more fish than the LME-mean), so the
   0-D q does **not** transfer; q is re-fit inside the gridded model (drops ~1–2 orders of magnitude).
4. **Time-varying, per-spectrum (U/V) fished-size selectivity** from the `_uv` parquet; the gravity/
   effort split uses the fishable biomass within each window.

## Scripts → ISIMIP3a pipeline stages

| Script | Role | Relates to stage |
|---|---|---|
| `07b_gridded_calibration_inputs.R` | **shim**: map stage-07 `best-fishing-parameters_*.parquet` (fmort_u/v → q_pel/q_ben) + the stage-04 U/V catch split into the per-region rds the gridded scripts read | bridges `07_estimating_best_vals` → gridded |
| `build_percell_bw.R` | per-cell, per-timestep biomass-weighted intercept/slope/T_exp/tob/export from the gridded GFDL-MOM6-COBALT2 netcdfs (+ gridded spin-up, 6× cycle of 1961-1980) | `01`, `03`, `07_setup_gridded` + `integrating_phyto(weighting="biomass")` |
| `build_center.R` | LME-center per-cell series on the biomass-weighted aggregate (preserves 0-D q) | new (post-`03`) |
| `gridded_calib.R` | re-fit (q_pel,q_ben) inside the gridded model (spin cached, BOBYQA on the transient aggregate) | new (post-`07_estimating_best_vals`) |
| `gridded_run.R` | the gridded dbpmr run: unfished spin → fished 1841-2010, annual gravity, warm-restart | dbpmr analog of `08_run_dbpm_gridded.py` |
| `run_gridded.sh` | resumable, retry-looped multi-region runner (RAM-disk TMPDIR, `caffeinate`) | orchestration |
| `plot_gridded.R` | per-region catch obs-vs-gridded + circumpolar maps | synthesis (`09`) |

## Southern Ocean validation (FAO 48/58/88 + LME 61 Antarctica)

| region | gridded q_pel | vs 0-D | corr (log10 catch) | level | verdict |
|---|---|---|---|---|---|
| FAO 58 Indian | 4.3e-5 | ÷83 | **+0.92** | 1.3× | good |
| FAO 48 Atlantic | 2.1e-4 | ÷15 | **+0.90** | 1.8× | good |
| FAO 88 Pacific | 2.5e-5 | ÷4 | −0.24 | ~0× | krill/ice-edge mismatch |
| LME 61 Antarctica | 4.5e-6 | ÷100 | +0.33 | ~0× | weak (krill/ice) |

**Finding:** the gridded DBPM reproduces conventional shelf/slope fisheries (Indian, Atlantic sectors:
r≈0.9, level within ~2×) but breaks down where the fishery is krill/ice-edge driven (Pacific sector,
coastal Antarctica) — the catch there tracks marginal-ice-zone dynamics outside the biomass-gravity
mechanics. A reportable model–data boundary, not a calibration artifact.

## Running

Requires: `dbpmr` installed (set `DBPMR_LIB`), the gridded GFDL netcdfs (`phyc/phypico/thetao/tob/
expc-bot/intpp`, 60 arcmin, in `gridded_nc/`), the `_uv` region parquets (`INPUT_PARQUET_DIR` or
`./dbpm_inputs_uv`), the 1° FAO-LME mask, and the 0-D calibration rds (`calib_A3/lme<L>.rds`).

```sh
# bridge the pipeline's 0-D calibration (stage 07) into the rds the gridded scripts read:
Rscript 07b_gridded_calibration_inputs.R --regions=158,148,.. --results=<07_out> --catch=<uv_split> --out=calib_A3
Rscript build_percell_bw.R <L..> --par=4          # per-cell biomass-weighted spatiotemporal inputs
Rscript build_center.R      <L..> --par=4          # LME-center -> percell_c_lme<L>.parquet
sh run_gridded.sh                                  # per region: gridded_calib.R (q refit) -> gridded_run.R (full grid)
Rscript plot_gridded.R                             # figures
```

`<CALIB_DIR>/lme<L>.rds` schema (produced by `07b`, read by the gridded scripts): `q_pel, q_ben`
(0-D catchabilities), `year, obs_pel, obs_ben` (observed U/V catch density — needed for the gridded
q refit), `region, corr_pel, corr_ben`.

Compute notes (laptop): dbpmr's warm-restart is file-I/O heavy — use a RAM-disk `TMPDIR` and (on
managed machines) exclude it from real-time AV scanning; `run_gridded.sh` is resumable at the region
level with a per-region retry loop.
