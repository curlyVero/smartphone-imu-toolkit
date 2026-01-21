# Smartphone IMU Toolkit (MATLAB + Python batch)

A lightweight toolkit for **smartphone IMU preprocessing and analysis**.  
This repository provides:

- **MATLAB tools** for interactive preprocessing, plotting, noise analysis, mounting-angle analysis, and KML generation.
- **Python batch scripts** that replicate the core logic of selected MATLAB tools for large-scale processing.

> **Equivalence**
> - `estimate_mounting_batch.py` implements the same mounting-angle analysis workflow as `mounting_Analysis.m`.
> - `process_imu.py` implements the same IMU preprocessing workflow as `process_ori_imu.m`.

---

## Contents

### MATLAB tools (interactive / analysis)
- `make_imuData.m`  
  Convert raw IMU logs into a consistent MATLAB-friendly format.
- `plot_imu_from_Uncal.m`  
  Plot IMU signals from *UncalAccel/UncalGyro* style logs.
- `plot_imu_timeline.m`  
  Plot IMU time series with a timeline layout for quick inspection.
- `process_ori_imu.m`  
  IMU preprocessing (e.g., basic cleaning, reformatting, time handling, and preparation for subsequent analysis).
- `imu_noise_analysis.m`  
  IMU noise analysis utilities (e.g., static segment statistics / Allan-style analysis depending on implementation).
- `mounting_Analysis.m`  
  Mounting-angle analysis (roll/pitch/yaw estimation from driving segments, correlation-based selection, etc.).
- `main_generate_kml.m`, `LineKml.m`  
  Export trajectories / segments to KML for visualization in Google Earth.

### Python tools (batch processing)
- `process_imu.py`  
  Batch preprocessing for many IMU files (same logic as `process_ori_imu.m`).
- `estimate_mounting_batch.py`  
  Batch mounting-angle estimation (same logic as `mounting_Analysis.m`).

---

## Repository structure (suggested)

```text
.
├─ matlab/
│  ├─ make_imuData.m
│  ├─ plot_imu_from_Uncal.m
│  ├─ plot_imu_timeline.m
│  ├─ process_ori_imu.m
│  ├─ imu_noise_analysis.m
│  ├─ mounting_Analysis.m
│  ├─ main_generate_kml.m
│  └─ LineKml.m
├─ python/
│  ├─ process_imu.py
│  └─ estimate_mounting_batch.py
└─ README.md
```

> If your repo currently places all scripts in the root directory, you can still use this README as-is.
> The folder layout above is recommended but not mandatory.

---

## Requirements

### MATLAB
- MATLAB R2018a+ is usually sufficient for most scripts.
- Toolboxes are generally not required unless your implementation calls specific toolbox functions.

### Python (batch scripts)
- Python 3.8+
- Typical dependencies (depending on your script implementation):
  - `numpy`
  - `pandas`
  - `scipy` (only if optimization / filtering is used)

Install dependencies (recommended):
```bash
pip install -r requirements.txt
```

If you do not have `requirements.txt`, you can install the basics:
```bash
pip install numpy pandas scipy
```

---

## Quick start

### 1) MATLAB: run a single file workflow
Open MATLAB, add the repo (or `matlab/`) to your path:

```matlab
addpath(genpath(pwd));
```

Typical workflows:

**(A) Preprocess one IMU file**
```matlab
process_ori_imu('path/to/input_imu_file', 'path/to/output_dir');
```

**(B) Mounting-angle analysis**
```matlab
mounting_Analysis('path/to/preprocessed_imu_file_or_dir');
```

**(C) Noise analysis**
```matlab
imu_noise_analysis('path/to/preprocessed_or_static_data');
```

> The exact function signatures may differ depending on your implementation.
> Please check the header comments at the top of each `.m` file for the authoritative usage.

---

### 2) Python: batch processing (recommended for many files)

#### (A) Batch IMU preprocessing — `process_imu.py`
This script is designed for processing **many files/folders** in one run.  
It follows the same preprocessing logic as `process_ori_imu.m`, but runs in batch mode.

Example:
```bash
python python/process_imu.py --input ./imu_ori --output ./imu_prepro
```

Common options (example, adjust to your script):
- `--input`: root directory of raw IMU files
- `--output`: output directory for preprocessed files
- `--pattern`: optional glob pattern for file matching
- `--recursive`: process recursively

#### (B) Batch mounting estimation — `estimate_mounting_batch.py`
This script is the batch version of `mounting_Analysis.m`.

Example:
```bash
python python/estimate_mounting_batch.py
```

By default, it typically uses a configuration block inside the script, for example:
- `DATA_ROOT = './imu_prepro'`
- `OUTPUT_FILE = 'mounting_summary.csv'`

After running, you should obtain a CSV summary such as:
- `mounting_summary.csv`

> Please check the header section in `estimate_mounting_batch.py` for the exact configuration keys and outputs.

---

## Input data conventions (recommended)

Because smartphone IMU logs vary across vendors/apps, this toolkit assumes you have (or can convert to) consistent columns such as:

- timestamp (sec)
- accelerometer (m/s²)
- gyroscope (rad/s)

If your raw logs are Android Sensor logs (e.g., UncalAccel/UncalGyro), use:
- `make_imuData.m` / `plot_imu_from_Uncal.m` to convert and verify,
then run preprocessing and analysis.

> **Tip:** Always confirm **units** before analysis.  
> Mixing `g` and `m/s²`, or `deg/s` and `rad/s`, will lead to incorrect estimates.

---

## Outputs

Depending on which tool you run, outputs may include:
- preprocessed IMU files (cleaned / reformatted)
- diagnostic plots (time series, segments, etc.)
- mounting-angle summary CSV (batch mode)
- KML files for trajectory visualization

---


---

## Contact / Issues

If you have any Questions or Comments, Please contact lwlLiu@wdu.edu.cn, Thanks!
