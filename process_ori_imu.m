% =========================================================================
% process_ori_imu - Interpolate and causally preprocess smartphone IMU "Uncal" logs.
%
% DESCRIPTION:
%   This batch tool processes all IMU CSV files under the folder ./imu_old/.
%   Each file is expected to contain Android-style "UncalAccel" and
%   "UncalGyro" messages with UTC milliseconds. The tool:
%     1) Splits accelerometer and gyroscope streams.
%     2) Merges duplicate timestamps by averaging samples.
%     3) Estimates sampling rate and maximum gaps (for reporting only).
%     4) Interpolates both streams onto a uniform time grid over the
%        overlapping time span (default 50 Hz).
%     5) Writes the interpolated "raw" IMU to ./imu/ (same filename).
%     6) Applies vehicle-oriented pseudo real-time preprocessing (strictly
%        causal): despiking + causal moving average, and writes results to
%        ./imu_preprocess/ (same filename).
%     7) Optionally plots raw vs preprocessed comparisons.
%
% USAGE:
%   process_ori_imu();             % default 50 Hz, no plots
%   process_ori_imu(50, true);     % 50 Hz, plot comparisons per file
%
% INPUTS:
%   targetFreq : Target output frequency in Hz (default 50).
%   doPlot     : Logical flag to generate comparison plots (default false).
%
% OUTPUT FILE FORMAT:
%   Each output file is a comma-separated text file with columns:
%     t_ms, AccX, AccY, AccZ, GyroX, GyroY, GyroZ
%   where t_ms is UTC milliseconds.
%
% =========================================================================

function process_ori_imu(targetFreq, doPlot)
% Process all smartphone IMU files under .\imu_old (vehicle-mounted scenario):
% 1) Read UncalAccel / UncalGyro messages
% 2) Remove duplicate timestamps (average samples with the same timestamp)
% 3) Estimate sampling rate and maximum gap (printed for statistics only)
% 4) Interpolate on the overlapping time span to a uniform time grid (default 50 Hz)
% 5) Write the interpolated "raw" IMU to .\imu\ with the same filename
% 6) Vehicle-oriented pseudo real-time preprocessing (causal, using only current/past samples):
%    - Despike (clamp + neighbor jump detection)
%    - Causal moving average (low-pass)
%   Output to .\imu_preprocess\ with the same filename
% 7) (Optional) Plot raw vs preprocessed comparison
%
% Usage:
%   process_ori_imu();             % default 50 Hz, no plots
%   process_ori_imu(50, true);     % 50 Hz, plot comparison for each file

    clc;
    close all;
    if nargin < 1
        targetFreq = 50;  % default output rate (Hz)
    end
    if nargin < 2
        doPlot = false;   % no plots by default
    end
    dt_target = 1 / targetFreq;

    inDir   = fullfile(pwd, 'imu_old');
    outDir  = fullfile(pwd, 'imu');              % interpolated raw output
    preDir  = fullfile(pwd, 'imu_preprocess');   % preprocessed output

    if ~exist(inDir, 'dir')
        error('Input folder does not exist: %s', inDir);
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    if ~exist(preDir, 'dir'), mkdir(preDir); end

    files = dir(inDir);
    for k = 1:numel(files)
        if files(k).isdir
            continue; % skip subdirectories
        end

        inFile   = fullfile(inDir,  files(k).name);
        outFile  = fullfile(outDir, files(k).name);
        preFile  = fullfile(preDir, files(k).name);

        fprintf('\n============================\n');
        fprintf('Processing file: %s\n', inFile);

        try
            process_single_file(inFile, outFile, preFile, dt_target, targetFreq, doPlot);
        catch ME
            warning('Failed to process %s: %s', files(k).name, ME.message);
        end
    end

    fprintf('\nAll files have been processed.\n');
end

%% ==== Process a single file ====
function process_single_file(inFile, outFile, preFile, dt_target, targetFreq, doPlot)

    %-- read file --%
    opts = detectImportOptions(inFile, ...
        'FileType', 'text', ...
        'Delimiter', ',', ...
        'PreserveVariableNames', true);

    % Force key variable types
    idMsg = strcmp(opts.VariableNames, 'MessageType');
    if any(idMsg)
        opts.VariableTypes{idMsg} = 'string';
    end
    idUtc = strcmp(opts.VariableNames, 'utcTimeMillis');
    if any(idUtc)
        opts.VariableTypes{idUtc} = 'double';
    end

    T = readtable(inFile, opts);

    % Check required columns
    requiredVars = {'MessageType','utcTimeMillis','MeasurementX','MeasurementY','MeasurementZ'};
    for i = 1:numel(requiredVars)
        if ~ismember(requiredVars{i}, T.Properties.VariableNames)
            error('File %s is missing required column: %s', inFile, requiredVars{i});
        end
    end

    msgType = string(T.MessageType);
    isAcc  = (msgType == "UncalAccel");
    isGyro = (msgType == "UncalGyro");

    if ~any(isAcc)
        warning('  No UncalAccel data found. Skipping this file.');
        return;
    end
    if ~any(isGyro)
        warning('  No UncalGyro data found. Skipping this file.');
        return;
    end

    Ta = T(isAcc,  :);
    Tg = T(isGyro, :);

    % --- absolute time (s) and measurements ---
    ta_raw = double(Ta.utcTimeMillis) * 1e-3;   % absolute time (s)
    ax_raw = Ta.MeasurementX;
    ay_raw = Ta.MeasurementY;
    az_raw = Ta.MeasurementZ;

    tg_raw = double(Tg.utcTimeMillis) * 1e-3;
    gx_raw = Tg.MeasurementX;
    gy_raw = Tg.MeasurementY;
    gz_raw = Tg.MeasurementZ;

    % --- remove duplicate timestamps (average) ---
    [ta_abs, ax, ay, az] = merge_duplicate_time(ta_raw, ax_raw, ay_raw, az_raw);
    [tg_abs, gx, gy, gz] = merge_duplicate_time(tg_raw, gx_raw, gy_raw, gz_raw);

    % --- sampling rate and max gap ---
    if numel(ta_abs) > 1
        dt_acc      = diff(ta_abs);
        freq_acc    = 1 / median(dt_acc);
        max_gap_acc = max(dt_acc);
    else
        freq_acc    = NaN;
        max_gap_acc = NaN;
    end

    if numel(tg_abs) > 1
        dt_gyro      = diff(tg_abs);
        freq_gyro    = 1 / median(dt_gyro);
        max_gap_gyro = max(dt_gyro);
    else
        freq_gyro    = NaN;
        max_gap_gyro = NaN;
    end

    fprintf('  Acc  estimated rate ~ %.2f Hz, max gap = %.4f s\n', freq_acc,  max_gap_acc);
    fprintf('  Gyro estimated rate ~ %.2f Hz, max gap = %.4f s\n', freq_gyro, max_gap_gyro);

    if freq_acc < 0.8*targetFreq || freq_gyro < 0.8*targetFreq
        warning('  Sampling rate is much lower than target %.2f Hz. Interpolation densifies but adds no new information.', targetFreq);
    end

    % --- overlapping time range of Acc/Gyro (absolute time) ---
    tStart = max(ta_abs(1),  tg_abs(1));
    tEnd   = min(ta_abs(end), tg_abs(end));

    if tEnd <= tStart
        warning('  No overlap between Acc/Gyro time spans. Skipping this file.');
        return;
    end

    % Use floor to avoid a tiny overshoot at the end
    nStep = floor((tEnd - tStart) / dt_target);
    t_grid_abs = tStart + (0:nStep).' * dt_target;   % absolute time grid (s)

    % --- linear interpolation on the absolute time grid ---
    accX = interp1(ta_abs, ax, t_grid_abs, 'linear');
    accY = interp1(ta_abs, ay, t_grid_abs, 'linear');
    accZ = interp1(ta_abs, az, t_grid_abs, 'linear');

    gyroX = interp1(tg_abs, gx, t_grid_abs, 'linear');
    gyroY = interp1(tg_abs, gy, t_grid_abs, 'linear');
    gyroZ = interp1(tg_abs, gz, t_grid_abs, 'linear');

    % --- output matrix: t_ms, AccX/Y/Z, GyroX/Y/Z ---
    outMat = [t_grid_abs*1000, accX, accY, accZ, gyroX, gyroY, gyroZ];

    %========================
    % First output: interpolated "raw" data -> imu/
    %========================
    fid = fopen(outFile, 'w');
    if fid < 0
        error('Cannot open output file: %s', outFile);
    end
    header = sprintf('utcms(%.2fHz),AccX,AccY,AccZ,GyroX,GyroY,GyroZ', targetFreq);
    fprintf(fid, '%s\n', header);
    fmt = '%.3f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n';
    fprintf(fid, fmt, outMat.');
    fclose(fid);
    fprintf('  Wrote interpolated file: %s (%d rows, target %.2f Hz)\n', outFile, size(outMat, 1), targetFreq);

    %========================
    % Second output: vehicle-oriented pseudo real-time preprocessing -> imu_preprocess/
    %========================
    acc_raw  = outMat(:, 2:4);  % Ax,Ay,Az
    gyro_raw = outMat(:, 5:7);  % Gx,Gy,Gz

    % 1) Pseudo real-time despiking (clamp + neighbor jump detection)
    acc_desp  = despike_imu(acc_raw,  'acc');
    gyro_desp = despike_imu(gyro_raw, 'gyro');

    % 2) Pseudo real-time causal moving average (low-pass)
    Nwin_acc  = max(3, round(0.10 / dt_target));  % ~0.10 s window
    Nwin_gyro = max(3, round(0.06 / dt_target));  % ~0.06 s window

    acc_filt  = causal_movmean(acc_desp,  Nwin_acc);
    gyro_filt = causal_movmean(gyro_desp, Nwin_gyro);

    outMat_pre = [outMat(:,1), acc_filt, gyro_filt];

    fid2 = fopen(preFile, 'w');
    if fid2 < 0
        error('Cannot open preprocessing output file: %s', preFile);
    end
    header2 = sprintf('utcms(%.2fHz),AccX,AccY,AccZ,GyroX,GyroY,GyroZ', targetFreq);
    fprintf(fid2, '%s\n', header2);
    fprintf(fid2, fmt, outMat_pre.');
    fclose(fid2);

    fprintf('  Wrote preprocessed file: %s (%d rows, causal denoising)\n', preFile, size(outMat_pre, 1));

    %========================
    % Optional: plot comparison
    %========================
    if doPlot
        % Use relative time (s) for visualization
        t_rel = (outMat(:,1) - outMat(1,1)) / 1000.0;
        [~, name, ext] = fileparts(inFile);
        figname = sprintf('%s%s', name, ext);
        plot_imu_compare(t_rel, acc_raw, gyro_raw, acc_filt, gyro_filt, figname);
    end
end

%% ==== Helper: merge duplicate timestamps (average per timestamp) ====
function [t_unique, x_u, y_u, z_u] = merge_duplicate_time(t, x, y, z)
    % Sort by time
    [t_sorted, idx] = sort(t);
    x = x(idx);
    y = y(idx);
    z = z(idx);

    % unique() returns t_unique and ic (mapping each sample to its unique index)
    [t_unique, ~, ic] = unique(t_sorted);

    % For each unique timestamp, average the corresponding X/Y/Z samples
    x_u = accumarray(ic, x, [], @mean);
    y_u = accumarray(ic, y, [], @mean);
    z_u = accumarray(ic, z, [], @mean);
end

%% ==== Helper: pseudo real-time despiking (vehicle-tuned thresholds) ====
function X_out = despike_imu(X_in, typeStr)
% X_in: [N x 3], three axes
% typeStr: 'acc' or 'gyro'
% Strictly causal: uses only the current and previous samples.
% Empirical settings for vehicle-mounted smartphone IMU:
%   - Acc unit is typically m/s^2; in vehicles it rarely exceeds +/-20 m/s^2
%   - Gyro thresholds below assume rad/s. If your data is deg/s, increase thresholds accordingly.

    X_out = X_in;
    N = size(X_in, 1);
    if N <= 1
        return;
    end

    switch lower(typeStr)
        case 'acc'
            % Physical clamp (m/s^2)
            absMax = 30;    % ~3 g
            % Neighbor jump threshold (m/s^2); larger jumps are treated as spikes
            jumpTh = 15;
        case 'gyro'
            % If your gyro data is in deg/s:
            % absMax = 500;
            % jumpTh = 200;
            % If your gyro data is in rad/s (default here):
            absMax = 10; 
            jumpTh = 5;
        otherwise
            absMax = 1e6;
            jumpTh = 1e6;
    end

    % 1) Simple physical clamp
    X_out = max(min(X_out, absMax), -absMax);

    % 2) Neighbor jump detection (only i and i-1), strictly causal
    for i = 2:N
        diffv = X_out(i,:) - X_out(i-1,:);
        if any(abs(diffv) > jumpTh)
            % Spike detected: replace with the previous value
            X_out(i,:) = X_out(i-1,:);
        end
    end
end

%% ==== Helper: causal moving average ====
function Y = causal_movmean(X, Nwin)
% Causal moving average: for each time i,
% average samples within [i-Nwin+1, ..., i]. Near the beginning, the window is shortened automatically.
% This is pseudo real-time (no future samples are used).
%
% X: [N x d]
% Nwin: window length (positive integer)

    if nargin < 2 || Nwin <= 1
        Y = X;
        return;
    end

    [N, d] = size(X);
    Y = zeros(N, d);

    % Running sum
    runSum = zeros(1, d);

    for i = 1:N
        runSum = runSum + X(i,:);

        if i > Nwin
            % Remove the oldest sample outside the window
            runSum = runSum - X(i-Nwin, :);
            M = Nwin;
        else
            % Not enough samples yet: window length = i
            M = i;
        end

        Y(i,:) = runSum / M;
    end
end

%% ==== Helper: plot raw vs preprocessed comparison ====
function plot_imu_compare(t, acc_raw, gyro_raw, acc_filt, gyro_filt, titleStr)
% t: [N x 1] relative time (s)
% acc_*:  [N x 3]
% gyro_*: [N x 3]

    figure('Name', ['IMU Acc Compare: ' titleStr], 'NumberTitle', 'off');
    axLabel = {'X', 'Y', 'Z'};
    for i = 1:3
        subplot(3,1,i);
        plot(t, acc_raw(:,i), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8); hold on;
        plot(t, acc_filt(:,i), 'b', 'LineWidth', 1.0);
        grid on;
        ylabel(['Acc ' axLabel{i} ' (m/s^2)']);
        if i == 1
            title(['Accelerometer: ' strrep(titleStr, '_', '\_')]);
            legend('Raw (interpolated)', 'Preprocessed', 'Location','best');
        end
        if i == 3
            xlabel('Time (s)');
        end
    end

    figure('Name', ['IMU Gyro Compare: ' titleStr], 'NumberTitle', 'off');
    for i = 1:3
        subplot(3,1,i);
        plot(t, gyro_raw(:,i), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8); hold on;
        plot(t, gyro_filt(:,i), 'r', 'LineWidth', 1.0);
        grid on;
        ylabel(['Gyro ' axLabel{i}]);
        if i == 1
            title(['Gyroscope: ' strrep(titleStr, '_', '\_')]);
            legend('Raw (interpolated)', 'Preprocessed', 'Location','best');
        end
        if i == 3
            xlabel('Time (s)');
        end
    end
end