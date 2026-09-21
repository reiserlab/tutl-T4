# Manuscript Figures — Handoff

This directory contains the MATLAB code that produces the **main figure** and
**supplementary figure** of the Reiser-lab T4/T5 *tutl⁻* manuscript. A small
amount of post-processing was done in Adobe Illustrator (panel A schematic on
each figure, minor cosmetic alignment) — see §6 below.

---

## 1. Entry points

| Figure | Run this | Calls |
|---|---|---|
| Main | `scripts/generate_manuscript_fig_main.m` | `generate_manuscript_fig('main')` |
| Supp | `scripts/generate_manuscript_fig_supp.m` | `generate_manuscript_fig('supp')` |

The two wrapper scripts are 5 lines each — each just sets up the path and calls
`generate_manuscript_fig` with the right mode argument.  All the logic lives in
three sub-figure functions:

| Function | Args | Purpose |
|---|---|---|
| `generate_manuscript_fig(mode)`        | `mode ∈ {'main','supp'}` | Combiner: builds 18 × 16 cm canvas, calls the two sub-figure functions, copies axes/annotations into a unified figure, exports PDF + PNG. |
| `generate_manuscript_fig_ef(axis_mode)` | `axis_mode ∈ {'pd','ortho'}` | Top-half flash sub-figure (1×11 trace tiles + amplitude summaries). |
| `generate_manuscript_fig_ds(speed_dps)` | `speed_dps ∈ {28, 56}` | Bottom-half bar-sweep sub-figure (ring of traces + polar plots + DSI/AR/Vm box plots). |

`mode='main'` ⇒ `axis_mode='pd'` and `speed_dps=56`.
`mode='supp'` ⇒ `axis_mode='ortho'` and `speed_dps=28`.

Outputs land in `<data_root>/manuscript_figures/`:
- `fig_main_ABCDEF_v3_<timestamp>.pdf` / `.png`
- `fig_supp_ABCDEF_v3_<timestamp>.pdf` / `.png`

`data_root` defaults to the value of `DATA_ROOT` set in the two wrapper scripts (see
§2), and can be overridden via the `opts.data_root` argument of the underlying functions.

---

## 2. Data layout (where to point the pipeline)

The experimental data are kept on **shared lab storage** (the `ttl_1DRF`
dataset). Each user typically works from a local copy of that tree on their
own machine. The figure scripts read everything from a single root directory
that you point at via the `DATA_ROOT` variable in
`scripts/generate_manuscript_fig_main.m` and `_supp.m`:

```matlab
DATA_ROOT = '<path to your local copy of the ttl_1DRF tree>';   % <-- edit for your setup
```

The pipeline expects the following layout under `DATA_ROOT`:

```
<DATA_ROOT>/
├── population_results/
│   ├── batch_results.mat                              ← required: late-batch results
│   └── ring_of_traces_cache_gauss_999_<speed>dps_abs.mat
│                                                     ← auto-generated on first run
│                                                       (one cache file per speed)
│
├── pre-bar-flash/
│   └── population_results/
│       └── batch_results_pre_bf.mat                   ← required: early-batch results
│
├── pre-bar-flash/{control,ttl}/{ON,OFF}/<exp_folders>/
│                                                     ← raw early-batch experiments
│                                                       (only read on cache miss; can be
│                                                       absent if the cache file is present)
│
├── <YYYY_MM_DD_HH_MM>/                                ← raw late-batch experiment folders
│   └── currentExp.mat                                   (only read by batch_analyze_1DRF.m
│                                                       when (re)building batch_results.mat;
│                                                       can be absent for the figure pipeline)
│
└── manuscript_figures/                                ← output: PDFs/PNGs land here (auto-created)
```

The files you need on disk to regenerate the figures:
- `<DATA_ROOT>/population_results/batch_results.mat`
- `<DATA_ROOT>/pre-bar-flash/population_results/batch_results_pre_bf.mat`
- the raw recording folders (`<DATA_ROOT>/<session>/` and
  `<DATA_ROOT>/pre-bar-flash/{control,ttl}/{ON,OFF}/<session>/`)

The two `.mat` files are produced by `src/analysis/protocol2/batch_analyze_1DRF.m`
and `batch_analyze_pre_bar_flash.m`.  The raw folders are still required: the
bar-sweep figure reads the raw logs to build its ring-of-traces cache (first run
per speed) and to compute the membrane-potential panel on every run, and the
square-flash figure (`generate_supp_fig_rf.m`) parses the raw flash data.

The `batch_results.mat` must contain the M6-aligned fields:
`pd_flash_m6_aligned`, `ortho_flash_m6_aligned`, `centroid_m6_rounded`,
`ortho_centroid_m6_rounded`, `pd_flash_baselines`, `ortho_flash_baselines`,
plus `dsi_vector` and `max_v_aligned`.  Re-run `batch_analyze_1DRF.m` (which
needs the raw experiment folders) if a pre-existing `batch_results.mat`
lacks these.

| Source                                | Path under `DATA_ROOT`                                                | Used by                |
|---|---|---|
| Late-batch results (M6-aligned)       | `population_results/batch_results.mat`                                | all sub-figures        |
| Early-batch pre-bar-flash results     | `pre-bar-flash/population_results/batch_results_pre_bf.mat`           | `_ds` (cache miss)     |
| Direction LUT                         | `src/analysis/protocol2/bar_lut.mat` (in this repo)                   | `_ds`                  |
| Ring-of-traces alignment cache        | `population_results/ring_of_traces_cache_gauss_999_<speed>dps_abs.mat`| `_ds` (auto-built)     |

Cell counts (with the canonical dataset):
- **Flash panels (late dataset, 25 cells):** T4 ctrl = 5, T4 *tutl⁻* = 5, T5 ctrl = 7, T5 *tutl⁻* = 8
- **Bar-sweep panels (combined 48 cells):** 23 early + 25 late

---

## 3. Final panel layout

Both figures share the same skeleton (top row = flash, bottom row = bar-sweep
ring + summaries).

### Main figure (PD-axis flash + 56 dps bar sweeps)

| Panel | Content |
|---|---|
| A | Schematic (Illustrator only, not from MATLAB) |
| B | Schematic (Illustrator only) |
| C | T4 + T5 PD-axis flash trace tiles (11 positions, M6-aligned, PD on left, ND on right) |
| D | T4 + T5 PD-axis depolarization amplitude vs position + per-position Wilcoxon asterisks + FWHM brackets |
| E | T4 ring-of-traces + polar (control + *tutl⁻*) |
| F | T5 ring-of-traces + polar (control + *tutl⁻*) |
| G | DSI box plots + Aspect-Ratio box plots |
| H | Vm pre vs during stim (box plots) |

### Supplementary figure (orthogonal-axis flash + 28 dps bar sweeps)

| Panel | Content |
|---|---|
| A | T4 + T5 **orthogonal-axis** flash tiles |
| B | T4 + T5 orthogonal-axis amplitude summary + asterisks + FWHM brackets |
| C | T4 ring-of-traces + polar at 28 dps |
| D | T5 ring-of-traces + polar at 28 dps |
| E | DSI + Aspect Ratio box plots |
| F | Vm pre and during stim box plots |

---

## 4. Analysis pipeline — flash sub-figure (`generate_manuscript_fig_ef.m`)

### 4.1 Trace alignment (M6 centroid)
Each cell's 11 flash traces (one per spatial position) are **spatially**
aligned before averaging across cells: the 68%-area centroid of the cell's
amplitude-by-position profile (rounded to the nearest position, "M6") is moved
to the central row 6, and the traces are reindexed accordingly.  Alignment is
precomputed and stored in `batch_results.results(k).pd_flash_m6_aligned` and
`.ortho_flash_m6_aligned`, computed by `compute_m6_centroid` +
`reindex_to_peak` in `src/analysis/protocol2/`.  Positions are on a common
2-pixel (2.5 deg) grid; the three sessions recorded at 1-pixel spacing
(`results(k).flash_pitch_px == 1`) contribute every second stimulus and only
to the central grid rows (`reindex_to_peak(..., step = 2)`).

### 4.2 PD-on-left ordering (`axis_mode='pd'` only)
Code stores positions ND→PD (1…11).  The PD branch of the function flips the
column order for display, so PD ends up at column 1 (`-5`) and ND at column 11
(`+5`):

```matlab
c_dep = c_dep(:, 11:-1:1);   % PD now at column 1, ND at column 11
positions   = -5:5;          % PD at -5 (left), ND at +5 (right)
```

Ortho branch (`axis_mode='ortho'`) does **not** flip.

### 4.3 Per-trial amplitude windows (sample indices, 10 kHz)
```
STIM_ONSET   = 5001
STIM_OFFSET  = 5801             % 80 ms flash
DEP_WINDOW   = [5001, 7000]     % 200 ms post-onset
HYP_WINDOW   = [6001, end]      % starts 100 ms after onset
DEP_PCTILE   = 99.9             % robust max
HYP_PCTILE   =  0.1             % robust min, then min(., 0)
```

`extract_robust_amps` (local function) returns N×11 matrices of dep / hyp
amplitudes per cell × position.

### 4.4 Per-position statistics
For each of the 11 spatial positions, a two-sided **Wilcoxon rank-sum** is
computed between control and *tutl⁻* per-cell amplitudes, no multiple-comparison
correction. Asterisks: `*` p<0.05, `**` p<0.01, `***` p<0.001.

### 4.5 FWHM analysis (panel D / B brackets)
In `add_fwhm_bars` (local function):
1. Compute the **group-mean amplitude profile** across positions (mean across cells, NaN-omitting).
2. Find the peak; linearly interpolate the half-max crossings on each side; FWHM = right_x − left_x.
3. Compute a **per-cell FWHM** the same way for each cell; compare control vs *tutl⁻* with two-sided Wilcoxon rank-sum (`ranksum`); annotate the bracket with `*`/`**`/`***`.

Bracket positions are drawn at the bottom of the amplitude axes (control bracket
lower, *tutl⁻* bracket above it).

**Computed FWHM values (current figures):**

| | T4 ctrl | T4 *tutl⁻* | T5 ctrl | T5 *tutl⁻* |
|---|---|---|---|---|
| **Main (PD axis), group-mean FWHM (positions)** | 5.23 | 6.83 | 5.85 | 6.73 |
| **Main, per-cell median (n)** | 5.07 (5) | 6.61 (5) | 5.21 (7) | 6.54 (8) |
| **Main, p (Wilcoxon)** | 0.0952 (n.s.) | | **0.0059 (\*\*)** | |
| **Supp (ortho), group-mean FWHM** | 8.34 | 9.51 | 6.17 | 7.17 |
| **Supp, per-cell median (n)** | 6.30 (5) | 8.00 (5) | 5.87 (7) | 6.46 (8) |
| **Supp, p (Wilcoxon)** | 0.2063 (n.s.) | | 0.0721 (n.s.) | |

### 4.6 Pooled rank-sum across positions
`compute_pooled_ranksum` does a sliding pooled comparison over 3-position
windows, used to compute the p-values overlaid in panel D / B as small
horizontal asterisk markers.  Pool centers come back from this function.

### 4.7 Trace rendering
- 10 kHz traces are block-averaged downsampled by `DS_FACTOR = 10` → 1 kHz for the PDF.
- X axis cropped to start at sample 3701 (≈130 ms pre-stimulus).
- Default y limits: `[-5, 25]` mV (ΔmV mode); amplitude summary y: `[0, 25]`.
- Stimulus is drawn as a green horizontal bar (full-opacity, 6–8% of y-range height).
- Tile separators between adjacent positions are thick white lines so each tile has a visible gap.

### 4.8 Color scheme

| Cell type | ctrl line | *tutl⁻* line | ctrl fill | *tutl⁻* fill |
|---|---|---|---|---|
| T4 (ON) | black `[0 0 0]` | red `[1 0 0]` | `[0.80 0.80 0.80]` | `[1 0.70 0.70]` |
| T5 (OFF) | dark grey `[0.4 0.4 0.4]` | burgundy `[0.8 0.2 0.2]` | `[0.70 0.70 0.70]` | `[0.90 0.60 0.60]` |
| Stimulus marker | green `[0.2 0.7 0.2]` | | | |

Fill alpha = 0.35; SEM band alpha = 0.20.

---

## 5. Analysis pipeline — bar-sweep sub-figure (`generate_manuscript_fig_ds.m`)

### 5.1 Speed configuration

| | Main (`speed_dps=56`) | Supp (`speed_dps=28`) |
|---|---|---|
| `ROW_OFFSET_LATE` | 16 | 0 |
| `PRE_BF_SPEED` | `'medium'` | `'slow'` |
| `FUNC_PAIR` | `[5 6]` | `[3 4]` |
| `FWHM_SAMPLES` (Gaussian kernel) | 2500 | 5000 |
| Ring x-axis zoom + thicker traces | yes | no |
| Polar overlay padding factor | 1.10 | 1.15 |
| Scale-bar duration | 500 ms | 1000 ms |

### 5.2 Combined 48-cell dataset
Cells from the **late** batch (`batch_results.mat`) and the **early**
"pre-bar-flash" batch (`batch_results_pre_bf.mat`) are concatenated:

- Late batch: rows `plot_order + ROW_OFFSET_LATE` from `bar_data` (16 directions × 2 speeds layout).
- Early batch: parsed via `parse_bar_data_pre_bf(f_data, v_data, PRE_BF_SPEED)`, then polarity-corrected via `correct_off_polarity_swap(..., FUNC_PAIR)`, then taken with `row_offset = 0`.

`plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16]` reorders the 16
directions into the canonical CCW polar sequence used by `find_pd_from_lut`
(right-hand-rule fixed in `ortho_pos_order`).

**Note:** the per-cell *voltage* extraction (Panel H / F) calls
`parse_bar_data_pre_bf(f_data_v, v_data_v)` *without* the speed argument,
because that helper extracts voltage waveforms common to both speeds.  This
asymmetry is preserved exactly in `generate_manuscript_fig_ds.m`.

### 5.3 Temporal alignment (Gaussian-convolution centroid)
For the ring-of-traces, each direction's trace is aligned to the centroid of
its own Gaussian-smoothed response near a reference anchor:
- Smooth with a Gaussian kernel of FWHM `FWHM_SAMPLES`.
- Within `±SEARCH_HALF = 6000` samples (±600 ms) of the stimulus midpoint, find the sample whose smoothed value is the centroid of mass.
- Reindex so that centroid lands at `REF_SAMPLE = DISPLAY_HALF + 1 = 14001`.
- Display window is symmetric ±14000 samples around the centroid (~2.8 s total at 10 kHz, downsampled by 10).

Cached output in `ring_of_traces_cache_gauss_999_<speed>dps_abs.mat`.  **Delete
the cache if you change parsing, polarity-correction, or alignment logic** —
otherwise the cached results from the old code will be reused.

### 5.4 Polar / DSI / Aspect-ratio metric
When `SUBTRACT_BASELINE = false` (default; ring traces shown in absolute
voltage), the polar plots and DSI / AR are still computed from
**baseline-subtracted** amplitudes. Baseline is the pre-stim window mean
(`BL_START:BL_END` = samples 1000–9000); response is the 99.9th percentile
of the stim window after subtraction.

Per-direction response amplitudes feed:
- **Polar plot:** group-mean ± SEM polygon, ctrl + *tutl⁻* overlaid.
- **DSI** = (PD − ND) / (PD + ND), where PD is the cell's preferred direction from the LUT.
- **Aspect ratio** = (PD + ND) / (OD₁ + OD₂), where ODs are the two directions orthogonal to PD.

**Per-direction Wilcoxon rank-sum** (ctrl vs *tutl⁻*) is also computed at all
16 PD-aligned directions and printed to the console + saved to
`<data_root>/manuscript_figures/fig_ds_direction_stats_<speed>dps.txt`.
P-values are reported both raw and FDR-corrected (Benjamini-Hochberg, q=0.05).
Counts at q<0.05 in the current figures:

|                | T4 (ON) | T5 (OFF) |
|---|---|---|
| **Main, 56 dps** | 7/16 directions | 5/16 directions |
| **Supp, 28 dps** | 2/16 directions | 4/16 directions |

DSI / AR / Vm box plots use 4 groups (T4 ctrl, T4 *tutl⁻*, T5 ctrl, T5 *tutl⁻*);
significance brackets are scaled to fit and stacked without overlap. The
brackets, p-values, and dot positions are computed *inside* the local
plotting helpers (`draw_boxplot_panel`, `draw_ar_panel`) — they are
deterministic from the per-cell input vectors.

### 5.5 Reference line
A thin black `−65 mV` horizontal line is drawn full-width behind every ring
trace.

### 5.6 Polar overlay annotations
- "control" and "*tutl⁻*" group labels in the upper-left of each ring panel; *tutl⁻* in italics, color-matched to data line.
- Three radial reference lines from the polar origin: PD (top), OD (right), ND (bottom), with text labels.
- 16 outward-pointing grey arrows, one per trace position × ring panel, drawn via `annotation(fig, 'arrow', ...)` at the upper-right corner of each trace tile.

---

## 6. Combiner layout (`generate_manuscript_fig.m`)

```
FIG_W = 18 cm,  FIG_H = 16 cm
Top portion    (flash):   y ∈ [0.56, 0.98]   (CD_S = 0.42, CD_O = 0.56)
Bottom portion (ring):    y ∈ [0.01, 0.54]   (EFGH_S = 0.53, EFGH_O = 0.01)
```

Axes from each sub-figure are copied with `copyobj`, then their `Position`
property is rescaled in y. Annotations (text boxes, arrows, lines) are copied
into the figure's annotation layer with the same y rescaling. Source figures
are closed after copying.

Export uses `exportgraphics` with `'ContentType', 'vector'` for PDF and
300 dpi PNG.  Note: the PDF includes a creation timestamp in metadata, so
**byte-equality of two PDFs is not guaranteed** even when the underlying data
is identical. Compare the data, not the file bytes.

---

## 7. Reproducing the figures from scratch

There are two phases: build the batch `.mat` files (one-time, slow), then
generate the figures (fast).

### 7.1 First-time setup — build the batch results

Edit the two paths near the top of `scripts/build_batch_results.m` for your
machine:

```matlab
DATA_ROOT     = '/path/to/ttl_1DRF';                  % your local copy
CIRCSTAT_PATH = '/path/to/CircStat2012a';             % CircStat toolbox
```

Then in MATLAB, from the repo root:

```matlab
run('scripts/build_batch_results.m')
```

This runs the two batch pipelines in sequence and produces:
- `<DATA_ROOT>/population_results/batch_results.mat` (late dataset, 25 cells)
- `<DATA_ROOT>/pre-bar-flash/population_results/batch_results_pre_bf.mat` (early dataset, 23 cells)

The two pipelines together take a few minutes. CircStat is needed only at
this stage; both batch functions error out with a clear message if the
toolbox is not on the path.

If a colleague has already produced these two `.mat` files (e.g. via
shared lab storage), you can skip 7.1 — just drop the files into the paths
above and proceed to 7.2.

### 7.2 Generate the figures

Edit `DATA_ROOT` near the top of
`scripts/generate_manuscript_fig_main.m` and
`scripts/generate_manuscript_fig_supp.m` so that it points at the same
local copy as in 7.1 (see §2 for the expected layout). Then:

```matlab
addpath(genpath('src'));
run('scripts/generate_manuscript_fig_main.m')   % main figure
run('scripts/generate_manuscript_fig_supp.m')   % supplementary
```

If you'd rather not edit the wrappers, you can call the underlying
function directly with an explicit data root:

```matlab
addpath(genpath('src'));
generate_manuscript_fig('main', struct('data_root', '/path/to/ttl_1DRF'));
generate_manuscript_fig('supp', struct('data_root', '/path/to/ttl_1DRF'));
```

The figure pipeline itself does not require CircStat — only the build step
in 7.1 does.

First figure-pipeline run rebuilds the ring-of-traces caches (slow, several
minutes).
Subsequent runs are fast.

If `batch_results.mat` is missing the M6 fields, re-run
`src/analysis/protocol2/batch_analyze_1DRF.m` on the 1DRF dataset.

---

## 8. Illustrator post-processing (manual)

The MATLAB PDFs are opened in Adobe Illustrator and lightly touched up:
- **Panel A schematic** is original Illustrator artwork (not produced by MATLAB).
- Minor cosmetic alignment / spacing of panel letters and inter-panel gaps.
- Possible tweaks to font weight or kerning for matching the manuscript template.

No data or quantitative annotations (numbers, p-values, FWHM brackets,
asterisks) are edited in Illustrator. If a future regeneration produces
visibly different stats text or layout, re-import the new PDF and reapply only
the schematic + cosmetic adjustments — do not re-edit numerics.

---

## 9. Conventions / gotchas for future edits

- **Local functions** (`compute_fwhm_positions`, `extract_robust_amps`, `add_fwhm_bars`, etc.) are defined at the bottom of each function file. They are not callable from outside the file.
- **PD-on-left flip happens after data extraction** (`(:,11:-1:1)`). Anything indexing into amplitude matrices in PD mode assumes column 1 = PD, column 11 = ND.
- **Statistics:** all per-position and FWHM tests are two-sided Wilcoxon rank-sum, **no multiple-comparison correction**. State this in the figure legend.
- **Color of *tutl⁻* labels** must match the data line color (T4 = red, T5 = burgundy). Italics for the genotype.
- **Cache invalidation:** if the M6 alignment, FWHM kernel, or pre-bar-flash parsing logic changes, delete `ring_of_traces_cache_gauss_999_*dps_abs.mat` to force rebuild.
- **`scripts/` is a new top-level directory** introduced by this manuscript pipeline.

---

## 10. File index (just this manuscript pipeline)

```
scripts/
  build_batch_results.m                       ← one-time setup: runs both batches
  generate_manuscript_fig_main.m              ← user-facing wrapper
  generate_manuscript_fig_supp.m              ← user-facing wrapper
  generate_manuscript_fig.m                   ← combiner function
  generate_manuscript_fig_ef.m                ← flash sub-figure function
  generate_manuscript_fig_ds.m                ← bar-sweep sub-figure function
  MANUSCRIPT_FIGURES.md                       ← this file

src/analysis/protocol2/
  batch_analyze_1DRF.m                        ← extended with M6 + dsi_vector
  batch_analyze_pre_bar_flash.m               ← new: early-batch (summer 2025) processor
  compute_m6_centroid.m                       ← new: 68%-area centroid helper
  reindex_to_peak.m                           ← new: row-shift helper
  pipeline/parse_bar_data_pre_bf.m            ← new: early-batch parser

src/analysis/helper/
  correct_off_polarity_swap.m                 ← new: pre-Oct-15-2025 fix
  find_pd_from_lut.m                          ← extended with ortho_pos_order

<data_root>/
  population_results/batch_results.mat        ← required for all figures
  pre-bar-flash/population_results/batch_results_pre_bf.mat
  population_results/ring_of_traces_cache_gauss_999_*.mat
  manuscript_figures/                         ← PDFs/PNGs land here
```
