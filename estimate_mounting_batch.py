#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""estimate_mounting_batch.py

Batch mounting-angle estimation for vehicle-mounted smartphone IMU data.

This script is the Python counterpart of `mounting_Analysis.m` and is intended
for batch processing over a folder tree.

What it does
------------
1) Recursively searches `DATA_ROOT` (default: ./imu_prepro) for IMU files.
2) Estimates roll/pitch using static leveling (gravity alignment).
3) Estimates yaw using a turning correlation method with a fallback for
   straight-line acceleration segments.
4) Writes a CSV report (`OUTPUT_FILE`) with per-file results.

Usage
-----
From the repository root:

    python estimate_mounting_batch.py

Edit the constants in the CONFIGURATION block below to match your data layout.

Dependencies
------------
- numpy
- pandas
- scipy
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import os
import glob
import sys

# ================= CONFIGURATION =================
DATA_ROOT = './imu_prepro'
OUTPUT_FILE = 'mounting_summary.csv'

# Thresholds
STATIC_VAR_THRES = 0.25       # Std Dev threshold for static detection
ACC_MOVE_THRES = 0.2          # Min Acc magnitude to consider "moving"
GYRO_TURN_MIN = 0.02          # Min YawRate to calculate Correlation (Turning)
GYRO_STRAIGHT_TARGET = 0.05   # Desired straight line threshold for fallback logic
CORR_THRESHOLD = 0.5          # Threshold to switch between Correlation and Acc alignment
# ===============================================

# Consistent with mounting_Analysis.m
GRAVITY_MAG = 9.80665  # m/s^2


def dcm_from_euler_xyz(roll: float, pitch: float, yaw: float) -> np.ndarray:
    """Direction Cosine Matrix matching MATLAB angle2dcm(roll,pitch,yaw,'XYZ')."""
    cx, sx = np.cos(roll), np.sin(roll)
    cy, sy = np.cos(pitch), np.sin(pitch)
    cz, sz = np.cos(yaw), np.sin(yaw)

    Rx = np.array([[1.0, 0.0, 0.0], [0.0, cx, -sx], [0.0, sx, cx]], dtype=float)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]], dtype=float)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]], dtype=float)

    return Rz @ Ry @ Rx


def _wrap_deg(angle_deg: float) -> float:
    """Wrap to [-180, 180)."""
    return (angle_deg + 180.0) % 360.0 - 180.0


def load_imu_file(path: str):
    """Load IMU file."""
    def _read_try(kwargs):
        try:
            df0 = pd.read_csv(path, on_bad_lines='skip', **kwargs)
            if df0 is None or df0.empty: return None
            df0 = df0.loc[:, ~df0.columns.astype(str).str.match(r'^Unnamed')]
            return df0
        except Exception:
            return None

    tries = [
        dict(sep=None, engine='python'),
        dict(delim_whitespace=True),
        dict(sep=','),
        dict(sep='\t'),
    ]

    df = None
    for kw in tries:
        df = _read_try(kw)
        if df is not None and not df.empty: break

    if df is None or df.empty:
        for kw in tries:
            kw2 = dict(kw)
            kw2['skiprows'] = 1
            df = _read_try(kw2)
            if df is not None and not df.empty: break

    if df is None or df.empty:
        raise RuntimeError('Empty or unreadable file')

    df.columns = [str(c).strip() for c in df.columns]

    needed = {'AccX', 'AccY', 'AccZ', 'GyroX', 'GyroY', 'GyroZ'}
    if needed.issubset(set(df.columns)):
        acc_b = df[['AccX', 'AccY', 'AccZ']].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
        gyro_b = df[['GyroX', 'GyroY', 'GyroZ']].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
    else:
        df_num = df.apply(pd.to_numeric, errors='coerce').dropna(how='all')
        mat = df_num.to_numpy(dtype=float)
        if mat.shape[1] < 7:
            raise RuntimeError(f'Numeric table has only {mat.shape[1]} columns (<7), cannot parse')
        acc_b = mat[:, 1:4]
        gyro_b = mat[:, 4:7]

    mask_valid = np.isfinite(acc_b).all(axis=1) & np.isfinite(gyro_b).all(axis=1)
    acc_b = acc_b[mask_valid]
    gyro_b = gyro_b[mask_valid]

    if acc_b.shape[0] < 10:
        raise RuntimeError('Too few valid IMU samples after parsing')

    return acc_b, gyro_b, int(acc_b.shape[0])


def estimate_mounting_angles(acc_b: np.ndarray, gyro_b: np.ndarray):
    """Estimate mounting angles with fallback logic for Yaw."""
    
    # --- 1) Static detection ---
    acc_norm = np.linalg.norm(acc_b, axis=1)
    win_size = 50
    acc_std = (
        pd.Series(acc_norm)
        .rolling(window=win_size, center=True, min_periods=1)
        .std()
        .fillna(np.inf)
        .to_numpy(dtype=float)
    )
    static_mask = acc_std < STATIC_VAR_THRES

    if int(static_mask.sum()) < 10:
        static_mean = acc_b.mean(axis=0)
        used_static = False
    else:
        static_mean = acc_b[static_mask].mean(axis=0)
        used_static = True

    target_g = np.array([0.0, 0.0, -GRAVITY_MAG], dtype=float)

    # --- 2) Roll/Pitch optimization ---
    def cost_rp(x: np.ndarray) -> float:
        roll, pitch = float(x[0]), float(x[1])
        R = dcm_from_euler_xyz(roll, -pitch, 0.0)
        e = (R @ static_mean.reshape(3,)) - target_g
        return float(np.dot(e, e))

    phi_init = np.arctan2(static_mean[1], static_mean[2])
    theta_init = np.arctan2(-static_mean[0], np.sqrt(static_mean[1] ** 2 + static_mean[2] ** 2))
    x0 = np.array([phi_init, theta_init], dtype=float)

    res = minimize(
        cost_rp,
        x0,
        method='Nelder-Mead',
        options={'xatol': 1e-8, 'fatol': 1e-12, 'maxiter': 4000, 'disp': False},
    )

    best_roll = float(res.x[0])
    best_pitch = float(res.x[1])
    rp_residual = float(res.fun) # Accuracy metric for RP

    R_level = dcm_from_euler_xyz(best_roll, -best_pitch, 0.0)

    # --- 3) Rotate all data to level frame ---
    acc_level = (R_level @ acc_b.T).T
    gyro_level = (R_level @ gyro_b.T).T

    acc_dyn_level = acc_level - target_g
    yaw_rate_meas = gyro_level[:, 2]

    # --- 4) Yaw Estimation ---
    
    # A. Calculate Correlation Grid (Primary Method)
    acc_dyn_norm = np.linalg.norm(acc_dyn_level, axis=1)
    # Mask for turning
    corr_mask = (acc_dyn_norm > ACC_MOVE_THRES) & (np.abs(yaw_rate_meas) > GYRO_TURN_MIN)
    
    # Default values
    best_yaw_deg = 0.0
    final_corr = 0.0
    yaw_method = 'None'
    yaw_std_deg = np.nan # Accuracy metric for Yaw
    
    # Run Grid Search
    if int(corr_mask.sum()) > 20 and float(np.std(yaw_rate_meas[corr_mask])) > 1e-6:
        valid_dyn = acc_dyn_level[corr_mask]
        valid_rate = yaw_rate_meas[corr_mask]
        xs = valid_dyn[:, 0]
        ys = valid_dyn[:, 1]
        rate_std = float(np.std(valid_rate))

        yaw_grid = np.arange(-180.0, 180.0 + 1e-9, 0.1)
        corrs = np.zeros_like(yaw_grid, dtype=float)

        if float(np.std(xs)) > 1e-6 or float(np.std(ys)) > 1e-6:
             for i, yaw_d in enumerate(yaw_grid):
                psi = np.deg2rad(float(yaw_d))
                # Rotate dyn acc to check lateral component
                lat_acc = -xs * np.sin(psi) + ys * np.cos(psi)
                if float(np.std(lat_acc)) > 1e-6:
                    corrs[i] = float(np.corrcoef(lat_acc, valid_rate)[0, 1])
                else:
                    corrs[i] = 0.0
        
        idx = int(np.argmax(corrs))
        max_corr_val = float(corrs[idx])
        best_yaw_corr_deg = float(yaw_grid[idx])
    else:
        max_corr_val = -1.0
        best_yaw_corr_deg = 0.0

    # B. Decision Logic based on Correlation Threshold
    if max_corr_val > CORR_THRESHOLD:
        # ---> CASE 1: High Correlation, use Turning Method
        best_yaw_deg = best_yaw_corr_deg
        final_corr = max_corr_val
        yaw_method = 'Correlation'
        yaw_std_deg = np.nan # Not applicable for correlation peak
        status = 'Success'
    else:
        # ---> CASE 2: Low Correlation, use Straight Line Acceleration Method
        # Logic: Find segments with high Acc but low YawRate. 
        # Assume acceleration is forward (Body X). 
        # Measure angle of acceleration vector in Level Frame (Ax, Ay).
        # yaw = arctan2(Ay, Ax)
        
        straight_mask = (acc_dyn_norm > ACC_MOVE_THRES) & (np.abs(yaw_rate_meas) < GYRO_STRAIGHT_TARGET)
        
        if int(straight_mask.sum()) > 10:
            valid_acc = acc_dyn_level[straight_mask]
            # Calculate angle of the acceleration vector
            # If Yaw=0, Acc should be [A, 0]. If Yaw=psi, Acc vector is rotated.
            # tan(psi) = Ay / Ax  -> psi = arctan2(Ay, Ax)
            inst_yaw_rad = np.arctan2(valid_acc[:, 1], valid_acc[:, 0])
            
            # Calculate mean angle (handling periodicity)
            mean_sin = np.mean(np.sin(inst_yaw_rad))
            mean_cos = np.mean(np.cos(inst_yaw_rad))
            best_yaw_rad = np.arctan2(mean_sin, mean_cos)
            
            best_yaw_deg = float(np.degrees(best_yaw_rad))
            final_corr = max_corr_val # Keep the low correlation for record
            yaw_method = 'Acceleration (Straight)'
            status = 'Success (Fallback)'
            
            # Accuracy: Std Dev of the instantaneous angles
            # Calculate difference from mean, wrap to pi, compute std
            diffs = np.arctan2(np.sin(inst_yaw_rad - best_yaw_rad), np.cos(inst_yaw_rad - best_yaw_rad))
            yaw_std_deg = float(np.degrees(np.std(diffs)))
            
        else:
            best_yaw_deg = 0.0
            final_corr = max_corr_val
            yaw_method = 'Failed'
            status = 'No turning or straight data'
            yaw_std_deg = np.nan

    # Final Output Construction
    angles_deg = np.array([
        -np.degrees(best_roll),
        -np.degrees(best_pitch),
        best_yaw_deg,
    ], dtype=float)

    # Wrap to [-180, 180)
    angles_deg = np.array([_wrap_deg(float(a)) for a in angles_deg], dtype=float)

    info = {
        'status': status,
        'yaw_method': yaw_method,
        'final_correlation': final_corr,
        'static_points': int(static_mask.sum()),
        'used_static_subset': bool(used_static),
        'turn_points': int(corr_mask.sum()),
        'rp_opt_success': bool(res.success),
        'rp_opt_fun': float(res.fun), # Residual cost
        'rp_opt_nit': int(getattr(res, 'nit', -1)),
        'static_mean_ax': float(static_mean[0]),
        'static_mean_ay': float(static_mean[1]),
        'static_mean_az': float(static_mean[2]),
        'yaw_std_deg': yaw_std_deg, # Precision info for Yaw (if Acc method)
        'rp_residual': rp_residual  # Precision info for Roll/Pitch (Optimization Cost)
    }

    return angles_deg, info


def process_file(path: str):
    """Load one file and run the algorithm."""
    try:
        acc_b, gyro_b, npts = load_imu_file(path)
        angles_deg, info = estimate_mounting_angles(acc_b, gyro_b)
        info['n_points'] = int(npts)
        status = str(info.get('status', 'Success'))
        return angles_deg, status, info
    except Exception as e:
        return None, f"Error: {str(e)}", {'status': f"Error: {str(e)}"}


def main():
    search_patterns = [
        os.path.join(DATA_ROOT, '**', '*IMU*'),
        os.path.join(DATA_ROOT, '**', '*.txt'),
        os.path.join(DATA_ROOT, '**', '*.csv')
    ]
    files = set()
    for p in search_patterns:
        files.update(glob.glob(p, recursive=True))
    files = sorted(list(files))

    print(f"Found {len(files)} files.")
    print(f"Processing (Threshold: Corr < {CORR_THRESHOLD} -> Use Straight Line Acc)...")

    results = []

    for i, fpath in enumerate(files):
        sys.stdout.write(f"\r[{i+1}/{len(files)}] Processing: {os.path.basename(fpath)}")
        sys.stdout.flush()

        angles, status, info = process_file(fpath)

        row = {
            'file_name': os.path.basename(fpath),
            'full_path': fpath,
            'status': status,
        }
        
        # Helper to safely get value
        def get_val(k, default=np.nan):
            return info.get(k, default)

        if angles is not None and 'Success' in status:
            row.update({
                'roll': float(angles[0]),
                'pitch': float(angles[1]),
                'yaw': float(angles[2]),
                'method': get_val('yaw_method', 'Unknown'),
                'final_correlation': float(get_val('final_correlation')),
                
                # --- Accuracy / Diagnostic Info ---
                'yaw_std_deg': float(get_val('yaw_std_deg')), # Precision for Yaw (Acc method)
                'rp_residual': float(get_val('rp_residual')), # Precision for Roll/Pitch (Cost)
                
                'static_points': int(get_val('static_points', 0)),
                'turn_points': int(get_val('turn_points', 0)),
                'rp_opt_success': bool(get_val('rp_opt_success', False)),
                'n_points': int(get_val('n_points', 0)),
            })
        else:
            row.update({
                'final_correlation': float(get_val('final_correlation')),
                'static_points': int(get_val('static_points', 0)),
                'n_points': int(get_val('n_points', 0)),
            })

        results.append(row)

    print("\nProcessing complete.")

    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_FILE, index=False)

    success = df[df['status'].str.contains('Success', na=False)]
    if not success.empty:
        print("\n--- Final Estimates (Mean over Success) ---")
        print(success[['roll', 'pitch', 'yaw']].mean())
        print(f"\nMethod Distribution:\n{success['method'].value_counts()}")

    print(f"Results saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()