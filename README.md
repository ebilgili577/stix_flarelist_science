# STIX Flarelist Science

Machine learning pipeline for solar flare localization using STIX (Spectrometer/Telescope for Imaging X-rays) detector data.

## Overview

This project provides tools for:
- **Fetching and processing** STIX flare data from the ESA Solar Orbiter mission
- **Training neural networks** to predict flare locations from raw detector counts
- **Using synthetic data** for model pre-training (generated externally)

## Installation

```bash
# Create conda environment
conda env create -f environment.yml
conda activate stix

# Install the package 
pip install -e .
```

## CLI Commands

After installation, three commands are available:

### `stix-fetch`
Fetch and process flare data from the STIX Data Center.

```bash
stix-fetch
```

### `stix-train`
Train localization models directly.

```bash
stix-train --experiment my_exp --mode mixed --col-count 12 --epochs 1000
```

Options:
- `--experiment` - Experiment name (required)
- `--mode` - Training mode: `real`, `syn`, `mixed`, `finetune`, or `all`
- `--col-count` - Number of detector subcollimators to use (default: 12)
- `--hidden-dims` - Model architecture, comma-separated (default: `1024,512,256,128,64`)
- `--lr` - Learning rate (default: 0.001)
- `--epochs` - Max epochs (default: 1000)
- `--batch` - Batch size (default: 512)
- `--patience` - Early stopping patience (default: 50)

Run `stix-train --help` for all options.

### `stix-submit`
Interactive SLURM job submission for training on a cluster.

```bash
stix-submit
```

Prompts for all parameters and submits a job via `sbatch`.

## Project Structure

```
stix_flarelist_science/
├── setup.py                 # Package configuration
├── environment.yml          # Conda environment
├── slurm/                   # SLURM job scripts
│   └── send_job.sh
├── stix_train/              # Training package
│   ├── cli.py               # stix-train entry point
│   ├── submit.py            # stix-submit entry point
│   ├── config.py            # Configuration and constants
│   ├── data.py              # Data loading and preprocessing
│   ├── model.py             # Neural network architecture
│   └── training.py          # Training loops
├── generate_flarelist/      # Flare data pipeline
│   ├── cli.py               # stix-fetch entry point
│   ├── flarelist_generate.py
│   └── ...
├── simulator/              # STIX data simulator (external)
│   └── [Clone stix-data-generator here]
├── data/                    # Data outputs
│   ├── 6_final/            # Final flarelist CSVs
│   ├── stix_data/          # Correction data files (phase_corr.sav, pitch_slit.sav, uv.sav)
│   └── synthetic/          # Synthetic data (generated externally)
└── experiments/             # Training outputs 
```

## Synthetic Data

This repository uses a previous version of the original STIX data simulator. For generating synthetic training data, please use the original and up-to-date version available at:

**[stix-data-generator](https://github.com/i4Ds/stix-data-generator)**

### Setup Instructions

1. **Clone the simulator**: Clone the `stix-data-generator` repository into the `simulator/` folder:
   ```bash
   git clone https://github.com/i4Ds/stix-data-generator.git simulator
   ```

2. **Place correction data**: Copy the required correction data files into `data/stix_data/`:
   - `phase_corr.sav`
   - `pitch_slit.sav`
   - `uv.sav`

3. **Generate synthetic data**: Use the simulator to generate synthetic data and save the output `.npz` files in the `data/synthetic/` directory. The training pipeline expects files named `sim_{n_samples}_{fov_big}.npz` (e.g., `sim_40000_4400.npz`).

## Training Modes

| Mode | Description |
|------|-------------|
| `real` | Train on real STIX data only |
| `syn` | Train on synthetic data only |
| `mixed` | Train on combined real + synthetic data |
| `finetune` | Pre-train on synthetic, fine-tune on real |
| `all` | Run all four modes sequentially |

## License

BSD 3-Clause License. See [LICENSE](LICENSE) for details.

Based on [stix_flarelist_science](https://github.com/hayesla/stix_flarelist_science) © 2023 Laura Hayes.
