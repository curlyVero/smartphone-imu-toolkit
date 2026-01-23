function out = imu_noise_analysis(imuFile, varargin)
% =========================================================================
% imu_noise_analysis - Estimate IMU noise metrics using overlapping Allan deviation.
%
% DESCRIPTION:
%   This function reads an interpolated/preprocessed IMU CSV file and reports
%   common noise indicators derived from Allan deviation, including:
%     - VRW  (Velocity Random Walk) for accelerometers
%     - ARW  (Angle Random Walk) for gyroscopes
%     - Bias instability (from the Allan deviation minimum)
%
% USAGE:
%   out = imu_noise_analysis('imu_preprocess/example.20IMU');
%   out = imu_noise_analysis('imu_preprocess/example.20IMU', ...
%       'Delimiter', ',', 'MakePlot', true);
%
% Key Features:
%   1. Robust Static Detection (removes handling noise).
%   2. Supports two input formats: Synced numeric and Android Uncal CSV.
%   3. Automatic White Noise region detection (slope -0.5 search).
%   4. Correlation Time (T) estimation using ACF:
%        - raw ACF on detrended signal (mostly reflects white noise -> very small T)
%        - bias ACF on low-pass / moving-mean bias series (recommended for GM bias)
%   5. Theoretical GM curve overlay to validate model fit.
%
% Output Metrics (Aligned with Continuous-Time PSD):
%   Model: b_dot = -(1/T)*b + w,   E[w^2] = q_c
% =========================================================================

close all;

% -------------------- Input Parsing --------------------
if nargin < 1 || (isstring(imuFile) && strlength(imuFile)==0) || (ischar(imuFile) && isempty(imuFile))
    error('Usage: out = imu_noise_analysis(imuFile, ...)');
end

p = inputParser;
p.addRequired('imuFile', @(x) ischar(x) || isstring(x));

p.addParameter('Delimiter', ',', @(x) ischar(x) || isstring(x));
p.addParameter('MakePlot', true, @(x) islogical(x) || isnumeric(x));

% Static Detection Parameters
p.addParameter('StaticDetect', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('GyroThr', 0.03, @(x) isnumeric(x) && isscalar(x) && x>0);    % rad/s
p.addParameter('AccThr', 0.15, @(x) isnumeric(x) && isscalar(x) && x>0);     % m/s^2 (|a|-g0)
p.addParameter('WinSec', 2.0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('MinStaticSec', 30, @(x) isnumeric(x) && isscalar(x) && x>0);

% Resampling / Downsampling
p.addParameter('ResampleStepMs', 10, @(x) isnumeric(x) && isscalar(x) && x>0); % For Uncal format
p.addParameter('TargetFsHz', 50, @(x) isnumeric(x) && isscalar(x) && x>0);     % Target Hz
p.addParameter('MaxAllanPoints', 100, @(x) isnumeric(x) && isscalar(x) && x>=30);

% White Noise Fitting (Sliding window on log10 tau)
p.addParameter('WN_WindowDecade', 0.6, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('WN_MinPts', 8, @(x) isnumeric(x) && isscalar(x) && x>=5);
p.addParameter('WN_MinR2', 0.90, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);

% NEW: Bias ACF extraction window (seconds)
% - This is critical: T for GM bias should be estimated on slow-varying bias,
%   not on raw series dominated by white noise.
p.addParameter('BiasACFWinSec', 120, @(x) isnumeric(x) && isscalar(x) && x>0);

p.parse(imuFile, varargin{:});

IMU_FILE  = char(p.Results.imuFile);
DELIM     = char(p.Results.Delimiter);
MAKE_PLOT = logical(p.Results.MakePlot);

StaticDetect  = logical(p.Results.StaticDetect);
GyroThr       = p.Results.GyroThr;
AccThr        = p.Results.AccThr;
WinSec        = p.Results.WinSec;
MinStaticSec  = p.Results.MinStaticSec;

TargetFsHz     = p.Results.TargetFsHz;
MaxAllanPoints = p.Results.MaxAllanPoints;
WN_Params = struct('win', p.Results.WN_WindowDecade, ...
                   'minPts', p.Results.WN_MinPts, ...
                   'minR2', p.Results.WN_MinR2);

BiasACFWinSec = p.Results.BiasACFWinSec;

g0 = 9.80665;

% -------------------- 1. Read & Preprocess --------------------
[t_ms, acc, gyro, input_format] = read_imu_any(IMU_FILE, DELIM, p.Results.ResampleStepMs);

% Convert to seconds, relative time
t = double(t_ms(:));
t = (t - t(1)) * 1e-3;
dt0 = median(diff(t));
fs0 = 1/dt0;

fprintf('--------------------------------------------------\n');
fprintf('Processing: %s\n', IMU_FILE);
fprintf('Format: %s | Fs_raw: %.2f Hz | Dur: %.2f s\n', input_format, fs0, t(end));

% -------------------- 2. Static Segment Detection --------------------
idx_use = true(size(t));
static_seg = [1 numel(t)];

if StaticDetect
    [static_seg, idx_use] = detect_static_segment(t, acc, gyro, g0, GyroThr, AccThr, WinSec, MinStaticSec);
    fprintf('Static Segment: [%.1f, %.1f] s (%.1f%% used)\n', ...
        t(static_seg(1)), t(static_seg(2)), 100*mean(idx_use));
else
    fprintf('StaticDetect: OFF (Using full data)\n');
end

t = t(idx_use);
acc = acc(idx_use,:);
gyro = gyro(idx_use,:);

% Recompute dt after cropping
dt = median(diff(t));
fs = 1/dt;

% -------------------- 3. Downsampling --------------------
downsample_factor = 1;
if fs > TargetFsHz * 1.2
    downsample_factor = max(1, round(fs / TargetFsHz));
    [t, acc, gyro] = downsample_avg(t, acc, gyro, downsample_factor);
    dt = median(diff(t));
    fs = 1/dt;
    fprintf('Downsampling: 1/%d -> Fs_new: %.2f Hz, Samples: %d\n', downsample_factor, fs, numel(t));
end

% -------------------- 4. Demean (Bias Removal) --------------------
acc0  = acc  - mean(acc, 1, 'omitnan');
gyro0 = gyro - mean(gyro, 1, 'omitnan');

% -------------------- 5. Correlation Time Estimation (ACF) --------------------
% IMPORTANT:
% - raw ACF on acc0/gyro0 is dominated by white noise, tends to yield tiny T (few samples).
% - bias ACF on low-pass / moving-mean bias series is recommended for GM bias model.
fprintf('Estimating Correlation Time using ACF ...\n');

% (5.1) Raw ACF (reference only)
T_acc_raw = zeros(3,1);
T_gyr_raw = zeros(3,1);
for k=1:3
    T_acc_raw(k) = estimate_corr_time_acf(acc0(:,k), fs);
    T_gyr_raw(k) = estimate_corr_time_acf(gyro0(:,k), fs);
end
T_acc_raw_avg = mean(T_acc_raw);
T_gyr_raw_avg = mean(T_gyr_raw);

% (5.2) Bias ACF (recommended)
% Extract "bias series" via moving mean
winN = max(3, round(BiasACFWinSec * fs));
acc_bias  = movmean(acc0,  winN, 1, 'omitnan');
gyro_bias = movmean(gyro0, winN, 1, 'omitnan');

T_acc_bias = zeros(3,1);
T_gyr_bias = zeros(3,1);
for k=1:3
    T_acc_bias(k) = estimate_corr_time_acf(acc_bias(:,k), fs);
    T_gyr_bias(k) = estimate_corr_time_acf(gyro_bias(:,k), fs);
end
T_acc_bias_avg = mean(T_acc_bias);
T_gyr_bias_avg = mean(T_gyr_bias);

fprintf('  -> Acc T_corr raw (avg):  %.4f s\n', T_acc_raw_avg);
fprintf('  -> Gyr T_corr raw (avg):  %.4f s\n', T_gyr_raw_avg);
fprintf('  -> Acc T_corr bias (avg): %.2f s  (BiasACFWinSec=%g s)\n', T_acc_bias_avg, BiasACFWinSec);
fprintf('  -> Gyr T_corr bias (avg): %.2f s  (BiasACFWinSec=%g s)\n', T_gyr_bias_avg, BiasACFWinSec);

% -------------------- 6. Allan Deviation (Overlapping) --------------------
fprintf('Computing Overlapping Allan Deviation ...\n');
[tauA, adevAx] = allan_adev_overlapping(acc0(:,1), dt, MaxAllanPoints);
[~,    adevAy] = allan_adev_overlapping(acc0(:,2), dt, MaxAllanPoints);
[~,    adevAz] = allan_adev_overlapping(acc0(:,3), dt, MaxAllanPoints);

[tauG, adevGx] = allan_adev_overlapping(gyro0(:,1), dt, MaxAllanPoints);
[~,    adevGy] = allan_adev_overlapping(gyro0(:,2), dt, MaxAllanPoints);
[~,    adevGz] = allan_adev_overlapping(gyro0(:,3), dt, MaxAllanPoints);

adevA_avg = (adevAx + adevAy + adevAz) / 3;
adevG_avg = (adevGx + adevGy + adevGz) / 3;

% -------------------- 7. Metrics Extraction --------------------
% (A) White Noise / Angle(Velocity) Random Walk
[fitA, idxWA] = fit_white_noise_window(tauA, adevA_avg, -0.5, WN_Params, 'ACC');
sigma_rw_acc = median(adevA_avg(idxWA) .* sqrt(tauA(idxWA)));  % m/s^2/sqrt(Hz)
VRW_avg = sigma_rw_acc * 60;                                  % m/s/sqrt(h)

if abs(fitA.slope - (-0.5)) > 0.2
    warning('[imu_noise_analysis] ACC: best window slope is far from -0.5. VRW may be unreliable. Consider longer / cleaner static data.');
end

[fitG, idxWG] = fit_white_noise_window(tauG, adevG_avg, -0.5, WN_Params, 'GYR');
sigma_rw_gyr = median(adevG_avg(idxWG) .* sqrt(tauG(idxWG)));  % rad/s/sqrt(Hz)
ARW_avg = (sigma_rw_gyr * 180/pi) * 60;                        % deg/sqrt(h)

if abs(fitG.slope - (-0.5)) > 0.2
    warning('[imu_noise_analysis] GYR: best window slope is far from -0.5. ARW may be unreliable. Consider longer / cleaner static data.');
end

% (B) Bias Instability (Sigma_bias)
[sigma_bi_acc, idxMinA] = bias_from_min(adevA_avg);
[sigma_bi_gyr, idxMinG] = bias_from_min(adevG_avg);

if idxMinG == numel(tauG)
    warning('[imu_noise_analysis] Gyroscope Allan deviation keeps decreasing; bias instability may be not observable with the current record length.');
end
if idxMinA == numel(tauA)
    warning('[imu_noise_analysis] Accelerometer Allan deviation keeps decreasing; bias instability may be not observable with the current record length.');
end

% -------------------- 8. Output Structure --------------------
out = struct();

% --- Primary Metrics (Average) ---
out.ACCEBIASSTD = sigma_bi_acc / 1e-5;                % mGal
out.GYROBIASSTD = (sigma_bi_gyr * 180/pi) * 3600;     % deg/h

out.VRW_mps_sqrth = VRW_avg;
out.ARW_deg_sqrth = ARW_avg;
out.AccNoiseDensity = sigma_rw_acc;
out.GyroNoiseDensity = sigma_rw_gyr;

% --- Correlation times (both raw and bias) ---
out.ACC_Tcorr_raw_s  = T_acc_raw_avg;
out.GYR_Tcorr_raw_s  = T_gyr_raw_avg;
out.ACC_Tcorr_bias_s = T_acc_bias_avg;
out.GYR_Tcorr_bias_s = T_gyr_bias_avg;

out.ACCELBTIME  = T_acc_bias_avg / 3600;              % hour (recommended for GM bias)
out.GYROBTIME   = T_gyr_bias_avg / 3600;              % hour (recommended for GM bias)

% --- GM Process Parameters (For Kalman Filter) ---
% Qc = 2 * sigma^2 / T, use bias T (recommended)
T_acc_use = max(T_acc_bias_avg, 1e-3);  % avoid division by zero
T_gyr_use = max(T_gyr_bias_avg, 1e-3);

out.ACC_BIAS_Qc = 2 * (sigma_bi_acc^2) / T_acc_use;
out.GYR_BIAS_Qc = 2 * (sigma_bi_gyr^2) / T_gyr_use;

% --- Stats ---
out.fs_Hz = fs;
out.static_duration_s = t(end) - t(1);
out.bias_acf_win_sec = BiasACFWinSec;

% -------------------- 9. Print Summary --------------------
fprintf('\n================ RESULTS ================\n');
fprintf('--- Accelerometer ---\n');
fprintf('VRW (White Noise) : %.3f m/s/sqrt(h)\n', VRW_avg);
fprintf('Bias Instability  : %.2f mGal\n', out.ACCEBIASSTD);
fprintf('T_corr raw (ACF)  : %.4f s\n', out.ACC_Tcorr_raw_s);
fprintf('T_corr bias (ACF) : %.2f s (%.4f h)  [BiasACFWinSec=%g s]\n', out.ACC_Tcorr_bias_s, out.ACCELBTIME, BiasACFWinSec);
fprintf('GM PSD (Qc)       : %.4e (m/s^2)^2/s   [use bias T]\n', out.ACC_BIAS_Qc);

fprintf('\n--- Gyroscope ---\n');
fprintf('ARW (White Noise) : %.3f deg/sqrt(h)\n', ARW_avg);
fprintf('Bias Instability  : %.2f deg/h\n', out.GYROBIASSTD);
fprintf('T_corr raw (ACF)  : %.4f s\n', out.GYR_Tcorr_raw_s);
fprintf('T_corr bias (ACF) : %.2f s (%.4f h)  [BiasACFWinSec=%g s]\n', out.GYR_Tcorr_bias_s, out.GYROBTIME, BiasACFWinSec);
fprintf('GM PSD (Qc)       : %.4e (rad/s)^2/s   [use bias T]\n', out.GYR_BIAS_Qc);
fprintf('=========================================\n');

% -------------------- 10. Plotting (With Theory Check) --------------------
if MAKE_PLOT
    plot_results(tauA, adevAx, adevAy, adevAz, adevA_avg, ...
                 sigma_rw_acc, sigma_bi_acc, T_acc_use, idxWA, ...
                 'Accelerometer', 'm/s^2');

    plot_results(tauG, adevGx, adevGy, adevGz, adevG_avg, ...
                 sigma_rw_gyr, sigma_bi_gyr, T_acc_use, idxWG, ...
                 'Gyroscope', 'rad/s');
end

end

% =========================================================================
% Helper: File IO
% =========================================================================
function [t_ms, acc, gyro, fmt] = read_imu_any(fname, delim, resampleStepMs)
    % Detect format by peeking at the header
    fid = fopen(fname,'r'); assert(fid>0,'Cannot open file: %s', fname);
    first = '';
    while true
        L = fgetl(fid);
        if ~ischar(L), break; end
        L = strtrim(L);
        if ~isempty(L), first = L; break; end
    end
    fclose(fid);

    is_uncal = contains(lower(first),'messagetype') || contains(lower(first),'uncal');

    if is_uncal
        fmt = 'uncal_csv';
        opts = detectImportOptions(fname,'FileType','text','Delimiter',delim);
        opts.ExtraColumnsRule = 'ignore';
        opts.EmptyLineRule = 'read';
        T = readtable(fname, opts);

        % Normalize column names
        v = lower(string(T.Properties.VariableNames));
        idx_msg = find(v=="messagetype",1);
        idx_t   = find(v=="utctimemillis",1);
        idx_mx  = find(v=="measurementx",1);

        if isempty(idx_mx), error('Uncal CSV format error: No MeasurementX found.'); end
        if isempty(idx_msg) || isempty(idx_t)
            error('Uncal CSV format error: Missing MessageType/utcTimeMillis.');
        end

        msg = string(T{:, idx_msg});
        tms = double(T{:, idx_t});
        raw = double(T{:, idx_mx:idx_mx+2});

        isAcc = contains(lower(msg),'accel');
        isGyr = contains(lower(msg),'gyro');

        [tA, aA] = make_unique_grid(tms(isAcc), raw(isAcc,:));
        [tG, gG] = make_unique_grid(tms(isGyr), raw(isGyr,:));

        fprintf('Uncal read: accel=%d, gyro=%d -> resample step=%g ms\n', numel(tA), numel(tG), resampleStepMs);

        % Align to uniform grid (overlapping time)
        t0 = max(min(tA), min(tG));
        t1 = min(max(tA), max(tG));
        grid_ms = (ceil(t0/resampleStepMs)*resampleStepMs : resampleStepMs : floor(t1/resampleStepMs)*resampleStepMs).';

        if isempty(grid_ms), error('No overlap between Acc and Gyro.'); end

        acc  = interp1(tA, aA, grid_ms, 'linear', 'extrap');
        gyro = interp1(tG, gG, grid_ms, 'linear', 'extrap');
        t_ms = int64(grid_ms);
    else
        fmt = 'synced_txt';
        try
            dat = readmatrix(fname); % Modern MATLAB
        catch
            dat = dlmread(fname, delim, 1, 0); % Legacy
        end
        if size(dat,2) < 7
            error('Synced format requires at least 7 columns: [t, ax, ay, az, gx, gy, gz]');
        end
        t_ms = int64(dat(:,1));
        acc  = dat(:,2:4);
        gyro = dat(:,5:7);
    end
end

function [tu, vu] = make_unique_grid(t, v)
    % Handle duplicate timestamps by averaging
    t = double(t(:));
    v = double(v);
    [t, idx] = sort(t);
    v = v(idx, :);

    ok = isfinite(t);
    t = t(ok);
    v = v(ok,:);

    [tu, ~, ic] = unique(t, 'stable');
    if numel(tu) < numel(t)
        vu = zeros(numel(tu), size(v,2));
        for k = 1:size(v,2)
            vu(:,k) = accumarray(ic, v(:,k), [], @mean);
        end
    else
        vu = v;
    end
end

% =========================================================================
% Helper: Signal Processing
% =========================================================================
function [T_corr] = estimate_corr_time_acf(x, fs)
    % Estimate correlation time T where ACF drops to 1/e
    % Robust to "immediate drop" cases dominated by white noise.

    x = x(:);
    x = x(isfinite(x));
    if numel(x) < 10
        T_corr = 0;
        return;
    end

    max_lag = round(fs * 3600 * 2); % up to 2 hours
    max_lag = min(max_lag, numel(x)-1);
    if max_lag < 2
        T_corr = 0;
        return;
    end

    [acf, lags] = xcorr(x, max_lag, 'coeff');

    pos_idx = lags >= 0;
    acf_p = acf(pos_idx);
    lags_p = lags(pos_idx);

    target = 1/exp(1);

    idx = find(acf_p <= target, 1, 'first');

    if isempty(idx)
        T_corr = lags_p(end) / fs;
        return;
    end

    if idx == 1
        % Already below threshold at lag=0? (can happen with numerical issues)
        T_corr = 0;
        return;
    end

    % Linear interpolation
    y1 = acf_p(idx-1); y2 = acf_p(idx);
    x1 = lags_p(idx-1); x2 = lags_p(idx);

    if abs(y2 - y1) < eps
        T_corr = x2 / fs;
    else
        lag_frac = x1 + (target - y1) * (x2 - x1) / (y2 - y1);
        T_corr = max(0, lag_frac / fs);
    end
end

function [seg, idx_use] = detect_static_segment(t, acc, gyro, g0, gyroThr, accThr, winSec, minStaticSec)
    N = numel(t);
    dt = median(diff(t));
    win = max(1, round(winSec/dt));

    gyro_norm = sqrt(sum(gyro.^2,2));
    acc_norm  = sqrt(sum(acc.^2,2));

    g_ref = median(acc_norm,'omitnan');

    is_static_sample = (gyro_norm < gyroThr) & (abs(acc_norm - g_ref) < accThr);

    k = ones(win,1);
    frac = conv(double(is_static_sample), k, 'same') / win;
    is_static_smooth = frac > 0.8;

    minLen = round(minStaticSec/dt);
    d = diff([false; is_static_smooth; false]);
    st = find(d==1);
    ed = find(d==-1)-1;
    lens = ed - st + 1;

    seg = [1 N];
    ok = false;

    if ~isempty(lens)
        good = find(lens >= minLen);
        if ~isempty(good)
            [~, im] = max(lens(good));
            i = good(im);
            seg = [st(i), ed(i)];
            ok = true;
        end
    end

    if ~ok
        warning('No static segment > %.1f s found. Using full record.', minStaticSec);
    end

    idx_use = false(N,1);
    idx_use(seg(1):seg(2)) = true;
end

function [t2, acc2, gyro2] = downsample_avg(t, acc, gyro, factor)
    N = numel(t);
    M = floor(N/factor);
    idx = reshape(1:(M*factor), factor, M);
    t2    = mean(t(idx),1).';
    acc2  = squeeze(mean(reshape(acc(1:M*factor,:), factor, M, 3), 1));
    gyro2 = squeeze(mean(reshape(gyro(1:M*factor,:), factor, M, 3), 1));
end

% =========================================================================
% Helper: Allan Deviation & Fitting
% =========================================================================
function [tau, adev] = allan_adev_overlapping(x, dt, nTau)
    x = x(:); x = x(~isnan(x));
    N = length(x);
    if N < 1000, error('Data too short for Allan analysis.'); end

    mMax = floor(N/3);
    mMin = 1;

    mList = unique(max(mMin, round(logspace(log10(mMin), log10(mMax), nTau))));

    tau  = zeros(size(mList));
    adev = zeros(size(mList));

    cs = [0; cumsum(x)];

    for i=1:length(mList)
        m = mList(i);
        if (N - 2*m) < 1
            tau(i)=NaN; 
            continue; 
        end

        avgs = (cs(1+m:N+1) - cs(1:N-m+1)) / m;

        if length(avgs) > m
            diffs = avgs(1+m:end) - avgs(1:end-m);
            adev(i) = sqrt(0.5 * mean(diffs.^2));
            tau(i) = m * dt;
        else
            tau(i) = NaN;
        end
    end

    valid = isfinite(tau) & isfinite(adev) & tau>0;
    tau = tau(valid);
    adev = adev(valid);
end

function [fit, idxBest] = fit_white_noise_window(tau, adev, targetSlope, params, tag)
    lt = log10(tau); 
    la = log10(adev);

    winDecade = params.win;
    Lmin = min(lt); Lmax = max(lt);
    step = 0.05;
    starts = Lmin : step : (Lmax - winDecade);

    bestScore = inf;
    fit = struct('slope',NaN,'R2',NaN);
    idxBest = 1:min(params.minPts, numel(tau));

    for s = starts
        e = s + winDecade;
        idx = find(lt>=s & lt<=e);
        if numel(idx) < params.minPts, continue; end

        X = lt(idx); Y = la(idx);
        P = polyfit(X, Y, 1);

        Yhat = polyval(P, X);
        SSres = sum((Y - Yhat).^2);
        SStot = sum((Y - mean(Y)).^2) + eps;
        R2 = 1 - SSres/SStot;

        score = abs(P(1) - targetSlope);
        if R2 < params.minR2, score = score + 10; end

        if score < bestScore
            bestScore = score;
            fit.slope = P(1);
            fit.R2 = R2;
            idxBest = idx;
        end
    end

    if isnan(fit.slope)
        idxBest = 1:min(10, numel(tau));
        fprintf('Warning: %s white noise fit failed, using start.\n', tag);
    end
end

function [sigma, idx] = bias_from_min(adev)
    [minVal, idx] = min(adev);
    sigma = minVal / 0.664;
end

% =========================================================================
% Helper: Plotting
% =========================================================================
function plot_results(tau, ax, ay, az, avg, sig_rw, sig_bi, T_corr, idxWN, name, unit)
    figure('Name', [name ' Allan Deviation'], 'Color', 'w');

    loglog(tau, ax, 'Color', [0.8 0.8 0.8]); hold on; grid on;
    loglog(tau, ay, 'Color', [0.8 0.8 0.8]);
    loglog(tau, az, 'Color', [0.8 0.8 0.8]);
    pAvg = loglog(tau, avg, 'b', 'LineWidth', 2);

    pWN = loglog(tau(idxWN), avg(idxWN), 'g.', 'MarkerSize', 15);

    ref_rw = sig_rw ./ sqrt(tau);
    ref_bi = sig_bi * ones(size(tau));

    loglog(tau, ref_rw, 'g--', 'LineWidth', 1.5);
    loglog(tau, ref_bi, 'b--', 'LineWidth', 1.5);

    xline(T_corr, 'k:', 'LineWidth', 1.5);

    xlabel('\tau (s)');
    ylabel(['Allan Deviation (' unit ')']);
    title([name ' Noise Analysis']);

    legend([pAvg, pWN], ...
       {'Average Axis', 'White Noise Fit'}, ...
       'Location', 'southwest');

    % ---------- Clear annotation (full metrics) ----------
    dim = [0.7 0.6 0.38 0.30];

    if contains(lower(name), 'acc')
        % ACC: sig_rw is noise density (m/s^2/sqrt(Hz)), VRW = sig_rw*60 (m/s/sqrt(h))
        VRW = sig_rw * 60;
        bias_mGal = sig_bi / 1e-5;                 % m/s^2 -> mGal
        Qc = 2 * (sig_bi^2) / max(T_corr, 1e-6);    % (m/s^2)^2/s
        str = {
            sprintf('VRW (White Noise): %.3f m/s/sqrt(h)', VRW)
            sprintf('Bias Instability : %.2f mGal', bias_mGal)
            sprintf('Corr Time        : %.2f s', T_corr)
            sprintf('GM PSD (Qc)      : %.4e (m/s^2)^2/s', Qc)
            };
    else
        % GYR: sig_rw is noise density (rad/s/sqrt(Hz)), ARW = sig_rw*180/pi*60 (deg/sqrt(h))
        ARW = (sig_rw * 180/pi) * 60;
        bias_degph = (sig_bi * 180/pi) * 3600;      % rad/s -> deg/h
        Qc = 2 * (sig_bi^2) / max(T_corr, 1e-6);    % (rad/s)^2/s
        str = {
            sprintf('ARW (White Noise): %.3f deg/sqrt(h)', ARW)
            sprintf('Bias Instability : %.2f deg/h', bias_degph)
            sprintf('Corr Time        : %.2f s', T_corr)
            sprintf('GM PSD (Qc)      : %.4e (rad/s)^2/s', Qc)
            };
    end

    annotation('textbox', dim, 'String', str, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'w', ...
        'EdgeColor', [0.7 0.7 0.7], 'FontName', 'Consolas', 'FontSize', 10);

end
