# Smartphone IMU Small Tools (MATLAB + Python)

A small collection of utilities for **smartphone IMU (Android UncalAccel/UncalGyro)** logs in vehicle scenarios.

This repo provides:
- MATLAB scripts for interactive inspection and single-file analysis
- Python scripts for large-scale batch processing

Features:
- merge raw accelerometer/gyroscope logs into a unified *Uncal-style* CSV
- visualize sampling instants and rate statistics
- interpolate to a uniform time grid and run **causal (pseudo real-time) preprocessing**
- estimate **IMU noise metrics** using overlapping Allan deviation
- estimate **mounting angles** (roll/pitch/yaw) for vehicle-mounted smartphones
- export a trajectory to **KML** for Google Earth

## Folder structure

Two folder conventions are supported.

### MATLAB workflow folders

Created automatically by `process_ori_imu.m`:

- `imu_old/`  
  Input folder. Place your original IMU CSV files here (Android Uncal format).

- `imu/`  
  Output folder. Interpolated “raw” IMU (uniform time grid).

- `imu_preprocess/`  
  Output folder. Causally preprocessed IMU (vehicle-oriented denoising).

### Python workflow folders

Used by `process_imu.py` (created automatically if missing):

- `imu_ori/`
  Input folder. Expected layout: `imu_ori/YYYY/DOY/*.IMU` (e.g., `2021/136/*.21IMU`).

- `imu_align/`
  Output folder. Interpolated “raw” IMU on a uniform time grid.

- `imu_prepro/`
  Output folder. Causally preprocessed IMU.

## Data formats

### 1) Uncal IMU CSV (Android style)

Expected columns:

```
MessageType,utcTimeMillis,MeasurementX,MeasurementY,MeasurementZ,BiasX,BiasY,BiasZ
```

- `MessageType`: `UncalAccel` or `UncalGyro`
- `utcTimeMillis`: UTC milliseconds (Unix epoch)
- `Measurement*`: accelerometer (m/s^2) or gyroscope (rad/s)
- `Bias*`: optional, can be zeros

`make_imuData.m` writes this format.

### 2) Interpolated IMU (written by `process_ori_imu.m`)

Columns:

```
t_ms, AccX, AccY, AccZ, GyroX, GyroY, GyroZ
```

## Quick start

### A. Merge separated ACC/GYR logs into Uncal CSV

If you have two separate files (accelerometer and gyroscope) in a simple numeric format:

```matlab
make_imuData('acc.txt', 'gyr.txt', 'imuData.csv');
plot_imu_from_Uncal('imuData.csv');
```

### B. Batch interpolate and preprocess all files under `imu_old/`

1) Put your Uncal CSV files into `./imu_old/`.

2) Run:

```matlab
process_ori_imu();          % default 50 Hz, no plots
process_ori_imu(50, true);  % 50 Hz, plot comparisons per file
```

Outputs will be written to:
- `./imu/` (interpolated raw)
- `./imu_preprocess/` (causal preprocessed)

### B2. Batch interpolate and preprocess with Python

Install Python dependencies:

```bash
pip install numpy pandas
```

Run batch processing (default expects `./imu_ori/YYYY/DOY/*.IMU`):

```bash
python process_imu.py --base_dir . --target_freq 50
```

Outputs will be written to:
- `./imu_align/YYYY/DOY/`
- `./imu_prepro/YYYY/DOY/`

Statistics and logs:
- `./imu_stats_detail.csv`
- `./imu_stats_summary.csv`
- `./imu_process.log`

### C. IMU noise analysis (Allan deviation)

Run Allan deviation on an interpolated/preprocessed file:

```matlab
out = imu_noise_analysis('imu_preprocess/example.20IMU', 'Delimiter', ',', 'MakePlot', true);
disp(out);
```

### D. Mounting angle estimation (vehicle scenario)

```matlab
result = mounting_Analysis('imu_preprocess/example.20IMU', 'Plot', true);
disp(result);
```

The returned angles are in degrees:
- `result.roll_deg`
- `result.pitch_deg`
- `result.yaw_deg`

### D2. Batch mounting angle estimation with Python

Install Python dependencies:

```bash
pip install numpy pandas scipy
```

By default, the script searches `./imu_prepro` recursively and writes
`mounting_summary.csv`:

```bash
python estimate_mounting_batch.py
```

If your preprocessed folder is different, update the `DATA_ROOT` constant at the
top of `estimate_mounting_batch.py`.

### E. Generate a KML trajectory

If your navigation/trajectory output file contains a header line with `GPSTime`, `Latitude`, and `Longitude` and then numeric lines:

```matlab
main_generate_kml('solution.txt', 'trajectory.kml', 'GNSS_INS_Trajectory');
```

Open the resulting `.kml` in Google Earth.

## Scripts overview

- `make_imuData.m`  
  Merge separate ACC/GYR logs into an Android-style Uncal CSV.

- `plot_imu_from_Uncal.m`  
  Inspect ACC/GYR sampling instants and basic rate statistics.

- `plot_imu_timeline.m`  
  Core plotting function used by `plot_imu_from_Uncal.m`.

- `process_ori_imu.m`  
  Batch interpolation to a uniform time grid + causal preprocessing for vehicle IMU.

- `imu_noise_analysis.m`  
  Overlapping Allan deviation + VRW/ARW/bias-instability reporting.

- `mounting_Analysis.m`  
  Roll/pitch from static leveling + yaw from correlation (grid search).

- `main_generate_kml.m`, `LineKml.m`  
  Export trajectory to KML (LineString).

- `process_imu.py` (Python)
  Batch interpolation to a uniform time grid + causal preprocessing.

- `estimate_mounting_batch.py` (Python)
  Batch mounting-angle estimation with correlation + straight-line fallback.

## Requirements

- MATLAB R2018b or later recommended.
- `detectImportOptions`, `readtable` are used for robust parsing.
- Optimization in `mounting_Analysis.m` uses `fminsearch` (built-in).

Python requirements:
- Python 3.8+ recommended
- `numpy`, `pandas` for `process_imu.py`
- `numpy`, `pandas`, `scipy` for `estimate_mounting_batch.py`

## Notes on units

- Accelerometer: **m/s^2**
- Gyroscope: **rad/s**  
  (`imu_noise_analysis.m` reports ARW in deg/sqrt(h) and bias instability in deg/h.)

## License

Choose a license before publishing (e.g., MIT, BSD-3-Clause, Apache-2.0).

