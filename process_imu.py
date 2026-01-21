#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""process_imu.py

Batch preprocessing for Android smartphone IMU logs (UncalAccel/UncalGyro).

What it does
------------
1) Recursively scans an input folder (default: ./imu_ori/YYYY/DOY/*IMU).
2) Splits each file into accelerometer and gyroscope streams by `MessageType`.
3) De-duplicates samples with identical timestamps (averaging duplicates).
4) Interpolates both streams onto a uniform time grid (default: 50 Hz).
5) Applies a simple vehicle-oriented causal preprocessing (spike suppression + EMA).
6) Writes two outputs per file:
   - aligned "raw" IMU on a uniform grid (default: ./imu_align/YYYY/DOY/)
   - causally preprocessed IMU (default: ./imu_prepro/YYYY/DOY/)
7) Exports per-file rate/gap statistics and a year-level summary.

Input format
------------
The input CSV is expected to contain at least these columns:

    MessageType, utcTimeMillis, MeasurementX, MeasurementY, MeasurementZ

`MessageType` should contain at least `UncalAccel` and `UncalGyro`.

Usage
-----
Typical usage (repo root):

    python process_imu.py --base_dir ./data --target_freq 50

This expects the following structure under --base_dir:

    imu_ori/YYYY/DOY/*.IMU

Outputs will be written to:

    imu_align/YYYY/DOY/
    imu_prepro/YYYY/DOY/

Dependencies
------------
- numpy
- pandas

Notes
-----
- The interpolation is done with NumPy only (no SciPy required).
- The preprocessing is designed to be causal/pseudo real-time (no look-ahead).
"""

from __future__ import annotations

import argparse
import logging
import os
from multiprocessing import Pool, cpu_count
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


LOGGER = logging.getLogger("process_imu")


def configure_logging(log_file: str) -> None:
    """Configure file-based logging.

    Note: Writing to the same log file from multiple processes is "best effort".
    For large deployments, consider redirecting stdout/stderr per process instead.
    """

    logging.basicConfig(
        filename=log_file,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )


def build_tasks(
    ori_root: str,
    align_root: str,
    prepro_root: str,
    target_freq: float,
) -> List[Tuple[str, str, str, float]]:
    """Build the multiprocessing task list.

    Each task is (in_file, out_align_dir, out_prepro_dir, target_freq).
    """

    tasks: List[Tuple[str, str, str, float]] = []

    if not os.path.isdir(ori_root):
        return tasks

    for year in sorted(os.listdir(ori_root)):
        year_dir = os.path.join(ori_root, year)
        if not os.path.isdir(year_dir):
            continue

        for doy in sorted(os.listdir(year_dir)):
            doy_dir = os.path.join(year_dir, doy)
            if not os.path.isdir(doy_dir):
                continue

            for fn in os.listdir(doy_dir):
                # Only process files whose extension ends with "IMU", e.g., ".20IMU"
                _, ext = os.path.splitext(fn)
                if not ext.endswith("IMU"):
                    continue

                in_file = os.path.join(doy_dir, fn)
                out_align_dir = os.path.join(align_root, year, doy)
                out_prepro_dir = os.path.join(prepro_root, year, doy)
                os.makedirs(out_align_dir, exist_ok=True)
                os.makedirs(out_prepro_dir, exist_ok=True)
                tasks.append((in_file, out_align_dir, out_prepro_dir, float(target_freq)))

    return tasks


def estimate_freq_and_gap(t_array: np.ndarray) -> Tuple[float, float]:
    """Estimate the sampling frequency and maximum gap based on timestamps (seconds)."""

    if len(t_array) <= 1:
        return float("nan"), float("nan")

    dt = np.diff(t_array)
    freq = 1.0 / float(np.median(dt))
    max_gap = float(np.max(dt))
    return float(freq), max_gap


def interp_1d_vec(t_src: np.ndarray, values: np.ndarray, t_tar: np.ndarray) -> np.ndarray:
    """Vectorized 1D linear interpolation for multi-channel signals.

    Parameters
    ----------
    t_src : (N,) strictly increasing source timestamps
    values : (N, d) source values
    t_tar : (M,) target timestamps

    Returns
    -------
    (M, d) interpolated values
    """

    values = np.asarray(values, dtype=float)
    t_src = np.asarray(t_src, dtype=float)
    t_tar = np.asarray(t_tar, dtype=float)

    d = int(values.shape[1])
    out = np.zeros((len(t_tar), d), dtype=float)
    for k in range(d):
        out[:, k] = np.interp(t_tar, t_src, values[:, k])
    return out


def preprocess_signal(X: np.ndarray, dt: float, mode: str) -> np.ndarray:
    """Causal preprocessing for vehicle-mounted smartphone IMU.

    This routine performs:
    1) Hard clipping (very conservative physical bounds)
    2) Slow EMA baseline estimation
    3) Spike detection (deviation from baseline OR large first difference)
    4) Fast EMA smoothing

    Parameters
    ----------
    X : (N, 3)
        Three-axis signal.
    dt : float
        Sampling interval in seconds.
    mode : {"acc", "gyro_rad"}
        Select parameter set.

    Returns
    -------
    (N, 3) filtered signal
    """

    if mode == "acc":
        # Units: m/s^2 (includes gravity projection). We mainly limit spikes and jerk.
        abs_max = 50.0
        spike_th = 8.0
        jerk_th = 80.0
        tau_base_s = 0.30
        tau_smooth_s = 0.12
    elif mode == "gyro_rad":
        # Units: rad/s
        abs_max = 6.0
        spike_th = 1.0
        jerk_th = 10.0
        tau_base_s = 0.40
        tau_smooth_s = 0.20
    else:
        raise ValueError(f"Unknown mode: {mode}")

    X = np.asarray(X, dtype=float)
    X = np.clip(X, -abs_max, abs_max)
    n, _ = X.shape

    # 1) Slow EMA baseline
    alpha_base = dt / (tau_base_s + dt)
    baseline = np.zeros_like(X)
    baseline[0] = X[0]
    for i in range(1, n):
        baseline[i] = baseline[i - 1] + alpha_base * (X[i] - baseline[i - 1])

    # 2) Spike suppression
    X_clean = X.copy()
    for i in range(1, n):
        dev = X[i] - baseline[i]
        jerk = (X[i] - X[i - 1]) / dt
        if np.any(np.abs(dev) > spike_th) or np.any(np.abs(jerk) > jerk_th):
            X_clean[i] = baseline[i]

    # 3) Fast EMA smoothing
    alpha_smooth = dt / (tau_smooth_s + dt)
    Y = np.zeros_like(X_clean)
    Y[0] = X_clean[0]
    for i in range(1, n):
        Y[i] = Y[i - 1] + alpha_smooth * (X_clean[i] - Y[i - 1])

    return Y


def save_csv(path: str, mat: np.ndarray, hz: float) -> None:
    """Save IMU matrix as CSV with a self-describing header."""

    header = f"utcms({hz:.2f}Hz),AccX,AccY,AccZ,GyroX,GyroY,GyroZ"
    np.savetxt(path, mat, delimiter=",", header=header, comments="", fmt="%.9f")


def process_single_file(
    in_file: str,
    align_dir: str,
    pre_dir: str,
    target_freq: float,
) -> Optional[Dict[str, object]]:
    """Process one input IMU CSV and write outputs.

    Returns a statistics dict (or None on failure).
    """

    fn = os.path.basename(in_file)
    LOGGER.info("Processing: %s", in_file)

    # Parse year and doy from path: .../imu_ori/YYYY/DOY/filename
    parts = in_file.split(os.sep)
    year = parts[-3] if len(parts) >= 3 else "unknown"
    doy = parts[-2] if len(parts) >= 2 else "unknown"

    try:
        df = pd.read_csv(in_file)
    except Exception as exc:
        LOGGER.error("%s: failed to read CSV: %s", fn, exc)
        return None

    required = [
        "MessageType",
        "utcTimeMillis",
        "MeasurementX",
        "MeasurementY",
        "MeasurementZ",
    ]
    for col in required:
        if col not in df.columns:
            LOGGER.error("%s: missing column '%s', skipped", fn, col)
            return None

    acc = df[df["MessageType"] == "UncalAccel"].copy()
    gyro = df[df["MessageType"] == "UncalGyro"].copy()
    if len(acc) == 0 or len(gyro) == 0:
        LOGGER.warning("%s: empty accel/gyro stream, skipped", fn)
        return None

    # Convert timestamps to seconds
    acc["t"] = acc["utcTimeMillis"].astype(float) * 1e-3
    gyro["t"] = gyro["utcTimeMillis"].astype(float) * 1e-3

    # De-duplicate: average samples with identical timestamps
    acc = (
        acc.groupby("t")[["MeasurementX", "MeasurementY", "MeasurementZ"]]
        .mean()
        .sort_index()
    )
    gyro = (
        gyro.groupby("t")[["MeasurementX", "MeasurementY", "MeasurementZ"]]
        .mean()
        .sort_index()
    )

    freq_acc, max_gap_acc = estimate_freq_and_gap(acc.index.values)
    freq_gyro, max_gap_gyro = estimate_freq_and_gap(gyro.index.values)
    LOGGER.info("%s: estimated rate Acc=%.3f Hz, Gyro=%.3f Hz", fn, freq_acc, freq_gyro)

    # Time overlap
    t_start = float(max(acc.index[0], gyro.index[0]))
    t_end = float(min(acc.index[-1], gyro.index[-1]))
    if t_end <= t_start:
        LOGGER.warning("%s: no time overlap, skipped", fn)
        return None

    dt = 1.0 / float(target_freq)
    n = int((t_end - t_start) / dt)
    if n < 2:
        LOGGER.warning("%s: too short overlap (n=%d), skipped", fn, n)
        return None

    t_grid = t_start + np.arange(n, dtype=float) * dt

    # Linear interpolation (NumPy only)
    acc_interp = interp_1d_vec(acc.index.values, acc.values, t_grid)
    gyro_interp = interp_1d_vec(gyro.index.values, gyro.values, t_grid)

    # Save aligned raw IMU
    align_mat = np.column_stack([t_grid * 1000.0, acc_interp, gyro_interp])
    out_align = os.path.join(align_dir, fn)
    save_csv(out_align, align_mat, target_freq)

    # Causal preprocessing
    acc_filt = preprocess_signal(acc_interp, dt, mode="acc")
    gyro_filt = preprocess_signal(gyro_interp, dt, mode="gyro_rad")
    pre_mat = np.column_stack([t_grid * 1000.0, acc_filt, gyro_filt])
    out_pre = os.path.join(pre_dir, fn)
    save_csv(out_pre, pre_mat, target_freq)

    return {
        "year": year,
        "doy": doy,
        "filename": fn,
        "freq_acc": float(freq_acc),
        "freq_gyro": float(freq_gyro),
        "max_gap_acc": float(max_gap_acc),
        "max_gap_gyro": float(max_gap_gyro),
        "n_acc": int(len(acc)),
        "n_gyro": int(len(gyro)),
    }


def _worker(args: Tuple[str, str, str, float]) -> Tuple[str, Optional[Dict[str, object]]]:
    """Multiprocessing worker wrapper."""

    in_file, out_align_dir, out_prepro_dir, target_freq = args
    stats = process_single_file(in_file, out_align_dir, out_prepro_dir, target_freq)
    return os.path.basename(in_file), stats


def write_stats_detail(path: str, rows: Iterable[Dict[str, object]]) -> None:
    """Write per-file statistics."""

    cols = [
        "year",
        "doy",
        "filename",
        "freq_acc",
        "freq_gyro",
        "max_gap_acc",
        "max_gap_gyro",
        "n_acc",
        "n_gyro",
    ]
    df = pd.DataFrame(list(rows), columns=cols)
    df.to_csv(path, index=False)


def generate_summary_report(detail_csv: str, summary_csv: str) -> None:
    """Generate year-level summary from the detail CSV."""

    if not os.path.exists(detail_csv):
        LOGGER.warning("Detail statistics not found: %s", detail_csv)
        return

    df = pd.read_csv(detail_csv)
    if df.empty:
        LOGGER.warning("Detail statistics is empty: %s", detail_csv)
        return

    summary = (
        df.groupby("year")
        .agg(
            n_files=("filename", "count"),
            avg_freq_acc=("freq_acc", "mean"),
            avg_freq_gyro=("freq_gyro", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(summary_csv, index=False)


def process_all_imu(
    ori_root: str,
    align_root: str,
    prepro_root: str,
    target_freq: float,
    stats_detail: str,
    stats_summary: str,
    workers: int,
) -> None:
    """Main entry point: build tasks, run multiprocessing, write statistics."""

    tasks = build_tasks(ori_root, align_root, prepro_root, target_freq)
    n_tasks = len(tasks)

    LOGGER.info("====================================")
    LOGGER.info("Found %d IMU files (extension *IMU)", n_tasks)
    LOGGER.info("CPU cores available: %d", cpu_count())
    LOGGER.info("Workers used: %d", workers)
    LOGGER.info("====================================")

    if not tasks:
        print("No IMU files found. Please check folder paths and file extensions.")
        return

    rows: List[Dict[str, object]] = []

    with Pool(processes=workers) as pool:
        for idx, (finished_name, stats) in enumerate(pool.imap_unordered(_worker, tasks), start=1):
            print(f"[{idx}/{n_tasks}] Done: {finished_name}", flush=True)
            if stats is not None:
                rows.append(stats)

    if rows:
        write_stats_detail(stats_detail, rows)
        generate_summary_report(stats_detail, stats_summary)

    LOGGER.info("All processing finished.")
    print("Finished. See log and statistics CSVs for details.")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Batch preprocess Android smartphone IMU logs (UncalAccel/UncalGyro)."
    )
    parser.add_argument(
        "--base_dir",
        type=str,
        default=".",
        help="Base directory that contains imu_ori/. Outputs are written under this base directory.",
    )
    parser.add_argument(
        "--ori_root",
        type=str,
        default=None,
        help="Input root. Default: <base_dir>/imu_ori",
    )
    parser.add_argument(
        "--align_root",
        type=str,
        default=None,
        help="Aligned output root. Default: <base_dir>/imu_align",
    )
    parser.add_argument(
        "--prepro_root",
        type=str,
        default=None,
        help="Preprocessed output root. Default: <base_dir>/imu_prepro",
    )
    parser.add_argument(
        "--target_freq",
        type=float,
        default=50.0,
        help="Target uniform sampling frequency in Hz.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, cpu_count() - 1),
        help="Number of worker processes (default: cpu_count()-1).",
    )
    parser.add_argument(
        "--stats_detail",
        type=str,
        default=None,
        help="Detail statistics CSV path. Default: <base_dir>/imu_stats_detail.csv",
    )
    parser.add_argument(
        "--stats_summary",
        type=str,
        default=None,
        help="Summary statistics CSV path. Default: <base_dir>/imu_stats_summary.csv",
    )
    parser.add_argument(
        "--log_file",
        type=str,
        default=None,
        help="Log file path. Default: <base_dir>/imu_process.log",
    )
    return parser.parse_args()


def main() -> None:
    """CLI entry."""

    args = parse_args()
    base_dir = os.path.abspath(args.base_dir)
    ori_root = args.ori_root or os.path.join(base_dir, "imu_ori")
    align_root = args.align_root or os.path.join(base_dir, "imu_align")
    prepro_root = args.prepro_root or os.path.join(base_dir, "imu_prepro")
    stats_detail = args.stats_detail or os.path.join(base_dir, "imu_stats_detail.csv")
    stats_summary = args.stats_summary or os.path.join(base_dir, "imu_stats_summary.csv")
    log_file = args.log_file or os.path.join(base_dir, "imu_process.log")

    os.makedirs(base_dir, exist_ok=True)
    configure_logging(log_file)

    process_all_imu(
        ori_root=ori_root,
        align_root=align_root,
        prepro_root=prepro_root,
        target_freq=float(args.target_freq),
        stats_detail=stats_detail,
        stats_summary=stats_summary,
        workers=int(args.workers),
    )


if __name__ == "__main__":
    main()
