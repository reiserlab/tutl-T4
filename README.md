# Nested RF and Direction Selectivity Protocols

Generates, runs and analyses two-stage receptive-field and direction-selectivity protocols
for patch-electrophysiology recordings of *Drosophila* T4/T5 neurons on the G4 LED arena.

**Documentation:** protocol and analysis pages are published at <https://reiserlab.github.io/tutl-T4/> (mirrored from `leburnett/reiser-documentation`).

**Version** README v1 · **Status** active · **Last verified** not yet verified

## What this does

- **Protocol 1** — presents pre-made flash grids to find the rough receptive-field location
  and preferred contrast of the recorded cell, giving a `peak_frame`.
- **Protocol 2** — generates a bespoke protocol centred on that `peak_frame`: high-resolution
  RF mapping with small flashing squares, plus direction selectivity from moving bars.
- **Analysis** — extracts direction-selectivity metrics and 2-D Gaussian receptive-field
  fits, and produces the manuscript figures.

**You will need:** a G4 LED arena, connected and calibrated ·
[`G4_Display_Tools`](https://github.com/leburnett/G4_Display_Tools) on the MATLAB path ·
MATLAB (Windows, for the rig) · recorded data — none is included in this repo.

## Quick start

Running an experiment needs the rig. To work with data you already have, start at
[Workflow](#workflow) step 3.

1. **Get the code**
   ```bash
   git clone https://github.com/leburnett/nested_RF_stimulus.git
   ```
2. **Add the source tree to the MATLAB path**
   ```matlab
   addpath(genpath('src'))
   ```
3. **Generate and run protocol 2**
   ```matlab
   generate_protocol2()
   ```
   Prompts for the `peak_frame` from protocol 1, the arena side, and the fly's age and
   strain; writes a timestamped experiment folder (`Patterns/`, `Functions/`,
   `currentExp.mat`) under `<matlabroot>\G4_Protocols\nested_RF_protocol2\`; runs it on the
   arena; prints the inverse peak frame for the opposite-contrast experiment.

## Repository map

| Directory | Contents |
|---|---|
| `src/protocol_generation/` | Builds and runs protocol 2 — [details](src/protocol_generation/README.md) |
| `src/stimulus_generation/` | Flash, bar and position-function generation for both protocols |
| `src/analysis/` | `process_protocol2.m` and the bar/flash/RF analyses |
| `src/preprocessing/` | Batch result generation and dashboard image export |
| `src/dashboard/` | Dash app for browsing per-cell results (pixi environment) |
| `scripts/` | Manuscript figures — [details](scripts/README.md) |
| `protocols/` | Pre-made protocol 1 `.g4p` files, by arena side and background level |
| `results/` | Pre-made patterns, position functions and generated outputs |
| `docs/` | [Protocol design and history](docs/protocol_background.md) |

## Workflow

The [Quarto guide](https://leburnett.github.io/reiser-documentation/Ephys/ephys_nested_rf.html)
is the canonical step-by-step version. In brief:

1. **Present protocol 1** — run a pre-made `.g4p` from `protocols/` (the folder matching the
   arena side and background level) via `G4_experiment_conductor`; read the `peak_frame` off
   the response plots.
2. **Generate and run protocol 2** — `generate_protocol2()`, as above.
3. **Process the data** — from *inside* the experiment folder:
   ```matlab
   cd('<path to>/nested_RF_protocol2/data/<timestamp>')
   process_protocol2()
   ```
   Writes `bar_results/` and `flash_results/`, plus figures.
   [Inputs and outputs](src/analysis/README.md).
4. **Browse and plot** — `cd src/dashboard && pixi run dashboard` (http://localhost:8051);
   manuscript figures from [`scripts/`](scripts/README.md).

## Key data types

| Output | Description | Units |
|---|---|---|
| `bar_results_*.mat` | Responses to 16 bar directions at 3 speeds, DSI, preferred direction | mV; deg |
| `rf_results_*.mat` | 2-D Gaussian RF fit, excitatory and inhibitory components | mV; arena px |
| Bar speeds | Sweeps presented at three speeds | 28, 56, 168 deg/s |

Stimulus grids and timing: [src/analysis/README.md](src/analysis/README.md).

## Conventions & gotchas

- **Bar stimuli are not regenerated.** Everything else is built de novo by
  `generate_protocol2()`, but the bar sweep and bar flash stimuli use pre-made patterns in
  `results/patterns/protocol2` — change their parameters by editing those patterns directly.
- `process_protocol2()` must be run from inside the experiment directory, not the repo.
- Protocol folders are named `yyyy_MM_dd_HH_mm` at generation time.
- Several analysis scripts still contain hardcoded data paths; edit before running elsewhere.

## Related

- [Full instructions and troubleshooting](https://leburnett.github.io/reiser-documentation/Ephys/ephys_nested_rf.html) — the canonical step-by-step guide.
- [`G4_Display_Tools`](https://github.com/leburnett/G4_Display_Tools) — arena control.

## Contact

Laura Burnett, Reiser Lab, HHMI Janelia Research Campus.
