% =========================================================================
% plot_imu_from_Uncal - Visualize timestamp alignment for UncalAccel/UncalGyro.
%
% DESCRIPTION:
%   Reads an "Uncal" IMU CSV file produced by make_imuData.m and plots the
%   sampling instants of accelerometer and gyroscope measurements to inspect
%   alignment and irregular sampling.
%
% USAGE:
%   plot_imu_from_Uncal('imuData.csv');
%
% INPUT:
%   imuFile: CSV with header
%     MessageType,utcTimeMillis,MeasurementX,MeasurementY,MeasurementZ,BiasX,BiasY,BiasZ
%
% OUTPUT:
%   A figure showing ACC and GYR sampling instants and basic rate statistics.
%
% DEPENDENCIES:
%   plot_imu_timeline.m
%
% =========================================================================

function plot_imu_from_Uncal(imuFile)
% Read imuData.txt (UncalAccel / UncalGyro) and plot sampling alignment
%
% imuData.txt format:
% MessageType,utcTimeMillis,MeasurementX,MeasurementY,MeasurementZ,BiasX,BiasY,BiasZ
% Example: plot_imu_from_Uncal("D:\WHU\25Sec\imu\demo\imu_old\p4l1136.20IMU")

    if nargin < 1
        error('Usage: plot_imu_from_Uncal(imuDataFile)');
    end

    % ---- read imuData ----
    T = readtable(imuFile, 'Delimiter', ',', 'FileType', 'text');

    % ---- separate ACC / GYR ----
    isAcc = strcmp(T.MessageType, 'UncalAccel');
    isGyr = strcmp(T.MessageType, 'UncalGyro');

    accUtcMs = T.utcTimeMillis(isAcc);
    gyrUtcMs = T.utcTimeMillis(isGyr);

    if isempty(accUtcMs) || isempty(gyrUtcMs)
        error('imuData.txt must contain both UncalAccel and UncalGyro');
    end

    % ---- sampling rate statistics ----
    rate_acc = imu_rate_stats(accUtcMs);
    rate_gyr = imu_rate_stats(gyrUtcMs);

    % ---- print info ----
    print_rate('ACC', rate_acc);
    print_rate('GYR', rate_gyr);

    % ---- plot ----
    plot_imu_timeline(accUtcMs, gyrUtcMs, rate_acc, rate_gyr);
end

% =========================================================
% Local utilities (kept here to avoid extra files)
% =========================================================
function rate = imu_rate_stats(utcMillis)
% Robust sampling rate estimation using median dt

    rate = struct('f', NaN, 'dt_min', NaN, 'dt_med', NaN, ...
                  'dt_max', NaN, 'n', numel(utcMillis));

    if numel(utcMillis) < 2
        return;
    end

    t = double(utcMillis(:)) / 1000.0;
    dt = diff(t);
    dt = dt(dt > 0);

    if isempty(dt)
        return;
    end

    rate.dt_min = min(dt);
    rate.dt_med = median(dt);
    rate.dt_max = max(dt);
    rate.f      = 1.0 / rate.dt_med;
end

function print_rate(tag, rate)
    if isnan(rate.f)
        fprintf('%s: insufficient or invalid timestamps (N=%d)\n', ...
            tag, rate.n);
        return;
    end

    fprintf('%s rate ~ %.3f Hz (N=%d, dt min/med/max = %.6f / %.6f / %.6f s)\n', ...
        tag, rate.f, rate.n, ...
        rate.dt_min, rate.dt_med, rate.dt_max);
end
