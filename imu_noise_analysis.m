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
%   The gyroscope unit is assumed to be rad/s in the input file.
%
% USAGE:
%   out = imu_noise_analysis('imu_preprocess/example.20IMU');
%   out = imu_noise_analysis('imu_preprocess/example.20IMU', ...
%       'Delimiter', ',', 'MakePlot', true, 'StaticFraction', 0.2);
%
% INPUT:
%   imuFile: Path to a CSV/text file with (at least) 7 numeric columns:
%       [t_ms, AccX, AccY, AccZ, GyroX, GyroY, GyroZ]
%     where t_ms is UTC milliseconds. This is the format written by
%     process_ori_imu.m in this repository.
%
% NAME-VALUE OPTIONS:
%   'Delimiter'      : Field delimiter in the file (default: ',').
%   'MakePlot'       : Whether to plot Allan deviation curves (default: true).
%   'StaticFraction' : Fraction of the earliest samples used as "static"
%                      data for Allan analysis (default: 0.2).
%
% OUTPUT:
%   out: Struct with per-axis estimates:
%     out.VRW_mps_sqrth      (1x3) accelerometer VRW in m/s/sqrt(h)
%     out.ACCEBIASSTD_mGal   (1x3) accelerometer bias instability in mGal
%     out.ARW_deg_sqrth      (1x3) gyroscope ARW in deg/sqrt(h)
%     out.GYROBIASSTD_degph  (1x3) gyroscope bias instability in deg/h
%
% NOTES:
%   - Allan deviation estimation requires sufficiently long static data.
%   - Bias instability is approximated by min(adev)/0.664.
%
% =========================================================================

    % ---- parse inputs ----
    if nargin < 1 || (isstring(imuFile) && strlength(imuFile) == 0) || (ischar(imuFile) && isempty(imuFile))
        error(['Usage: out = imu_noise_analysis(imuFile, ''Delimiter'', '','', ''MakePlot'', true)\n' ...
               'Example: out = imu_noise_analysis(''imu_preprocess/example.20IMU'');']);
    end

    p = inputParser;
    p.addRequired('imuFile', @(x) ischar(x) || isstring(x));
    p.addParameter('Delimiter', ',', @(x) ischar(x) || isstring(x));
    p.addParameter('MakePlot', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('StaticFraction', 0.2, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    p.parse(imuFile, varargin{:});

    IMU_FILE = char(p.Results.imuFile);
    DELIM    = char(p.Results.Delimiter);
    MAKE_PLOT = logical(p.Results.MakePlot);
    staticFrac = p.Results.StaticFraction;

    format long

    g0 = 9.80665;

    % -------- robust read (works for *.IMU written by process_ori_imu.m) --------
    opts = detectImportOptions(IMU_FILE, 'FileType','text', 'Delimiter', DELIM);
    opts.ExtraColumnsRule = 'ignore';
    opts.EmptyLineRule    = 'read';
    opts = setvartype(opts, 1:7, 'double');
    T = readtable(IMU_FILE, opts);

    t_ms = T{:,1};
    acc  = T{:,2:4};     % m/s^2
    gyro = T{:,5:7};     % rad/s

    t  = (t_ms - t_ms(1)) * 1e-3;   % s
    dt = median(diff(t));
    fs = 1/dt;

    % ---- use only an early portion of the data (assumed static) ----
    N = numel(t);
    N_use = max(10, floor(N * staticFrac));

    fprintf('File: %s\n', IMU_FILE);
    fprintf('fs = %.6f Hz, duration = %.2f s\n', fs, t(end));
    fprintf('Using first %.0f%% of samples as static: %d / %d (%.1f s)\n', ...
        staticFrac*100, N_use, N, t(N_use) - t(1));

    t    = t(1:N_use);
    acc  = acc(1:N_use, :);
    gyro = gyro(1:N_use, :);

    % Gyro in deg/s is used only for reporting (ARW in deg/sqrt(h), bias in deg/h)
    gyro_degps = gyro * 180/pi;

    % De-mean
    acc0      = acc        - mean(acc,1,'omitnan');
    gyro0_rad = gyro       - mean(gyro,1,'omitnan');

    % Allan deviation (overlapping)
    [tauA, adevAx] = allan_adev_overlapping(acc0(:,1), dt);
    [~,   adevAy]  = allan_adev_overlapping(acc0(:,2), dt);
    [~,   adevAz]  = allan_adev_overlapping(acc0(:,3), dt);
    adevA_avg = (adevAx + adevAy + adevAz)/3;

    [tauG, adevGx] = allan_adev_overlapping(gyro0_rad(:,1), dt);
    [~,   adevGy]  = allan_adev_overlapping(gyro0_rad(:,2), dt);
    [~,   adevGz]  = allan_adev_overlapping(gyro0_rad(:,3), dt);
    adevG_avg = (adevGx + adevGy + adevGz)/3;

    % Slope tolerance for selecting the -1/2 region (white noise)
    tolSlope = 0.12;

    % ---------- per-axis outputs ----------
    out = struct();
    out.VRW_mps_sqrth       = zeros(1,3);
    out.ACCEBIASSTD_mGal    = zeros(1,3);
    out.ARW_deg_sqrth       = zeros(1,3);
    out.GYROBIASSTD_degph   = zeros(1,3);

    adevA_all = {adevAx, adevAy, adevAz};
    adevG_all = {adevGx, adevGy, adevGz};

    for k = 1:3
        % Accelerometer white noise (-0.5 slope) -> sigma_rw_acc (m/s^2/sqrt(Hz))
        idxW = pick_slope_region(tauA, adevA_all{k}, -0.5, tolSlope, 'smallTau');
        sig_rw_acc = median(adevA_all{k}(idxW).*sqrt(tauA(idxW)));
        out.VRW_mps_sqrth(k) = sig_rw_acc * 60; % m/s/sqrt(h)

        % Accelerometer bias instability from Allan minimum
        sig_bi_acc = min(adevA_all{k}) / 0.664; % m/s^2
        out.ACCEBIASSTD_mGal(k) = sig_bi_acc / 1e-5; % mGal (1 mGal = 1e-5 m/s^2)

        % Gyro white noise (-0.5 slope) -> sigma_rw (rad/s/sqrt(Hz)) -> deg/s/sqrt(Hz)
        idxWg = pick_slope_region(tauG, adevG_all{k}, -0.5, tolSlope, 'smallTau');
        sig_rw_gyr_rad = median(adevG_all{k}(idxWg).*sqrt(tauG(idxWg)));
        sig_rw_gyr_deg = sig_rw_gyr_rad * 180/pi;
        out.ARW_deg_sqrth(k) = sig_rw_gyr_deg * 60; % deg/sqrt(h)

        % Gyro bias instability from Allan minimum
        sig_bi_gyr_rad = min(adevG_all{k}) / 0.664; % rad/s
        sig_bi_gyr_degps = sig_bi_gyr_rad * 180/pi;
        out.GYROBIASSTD_degph(k) = sig_bi_gyr_degps * 3600; % deg/h
    end

    % ---------- convenience printout based on the average curve ----------
    idxW_A = pick_slope_region(tauA, adevA_avg, -0.5, tolSlope, 'smallTau');
    sigma_rw_acc = median(adevA_avg(idxW_A).*sqrt(tauA(idxW_A)));       % m/s^2/sqrt(Hz)
    [sigma_bi_acc, idxMinA] = bias_from_min(adevA_avg);                % m/s^2
    tau_rw_A = median(tauA(idxW_A));
    tau_bi_A = tauA(idxMinA);

    idxW_G = pick_slope_region(tauG, adevG_avg, -0.5, tolSlope, 'smallTau');
    sigma_rw_gyr_rad = median(adevG_avg(idxW_G).*sqrt(tauG(idxW_G)));   % rad/s/sqrt(Hz)
    [sigma_bi_gyr_rad, idxMinG] = bias_from_min(adevG_avg);            % rad/s
    tau_rw_G = median(tauG(idxW_G));
    tau_bi_G = tauG(idxMinG);

    fprintf('\n=== Accelerometer (Average X-Y-Z) ===\n');
    fprintf('sigma_rw = %.6g m/s^2/sqrt(Hz) = %.2f m/s/sqrt(h)\n', sigma_rw_acc, sigma_rw_acc*60);
    fprintf('sigma_bi = %.6g m/s^2 = %.2f ug = %.2f mGal\n', sigma_bi_acc, sigma_bi_acc/g0*1e6, sigma_bi_acc/1e-5);

    fprintf('\n=== Gyroscope (Average X-Y-Z) ===\n');
    fprintf('sigma_rw = %.6g rad/s/sqrt(Hz) = %.2f deg/sqrt(h)\n', sigma_rw_gyr_rad, (sigma_rw_gyr_rad*180/pi)*60);
    fprintf('sigma_bi = %.6g rad/s = %.2f deg/h\n', sigma_bi_gyr_rad, (sigma_bi_gyr_rad*180/pi)*3600);

    % ---------- plots ----------
    if MAKE_PLOT
        plot_acc(tauA, adevAx, adevAy, adevAz, adevA_avg, sigma_rw_acc, sigma_bi_acc, tau_rw_A, tau_bi_A);
        plot_gyr(tauG, adevGx, adevGy, adevGz, adevG_avg, sigma_rw_gyr_rad, sigma_bi_gyr_rad, tau_rw_G, tau_bi_G);
    end
end


% ===================== plots =====================
function plot_acc(tau, ax, ay, az, aavg, sigma_rw, sigma_bi, tau_rw, tau_bi)
    figure('Name','Accelerometer Allan Deviation'); clf;
    loglog(tau, ax, 'LineWidth', 1.0); hold on; grid on;
    loglog(tau, ay, 'LineWidth', 1.0);
    loglog(tau, az, 'LineWidth', 1.0);
    loglog(tau, aavg,'LineWidth', 1.6);

    xlim([1e-2 1e4]);
    set(gca,'XTick',[1e-2 1e-1 1 10 100 1e3 1e4]);
    
    ref_rw  = sigma_rw ./ sqrt(tau);
    ref_bi  = sigma_bi * ones(size(tau));
    loglog(tau, ref_bi, '--', 'LineWidth', 1.2);
    loglog(tau, ref_rw,  '--', 'LineWidth', 1.2);

    y_rw = sigma_rw / sqrt(tau_rw);
    y_bi = sigma_bi;
    xline(tau_rw, '--'); yline(y_rw, '--');
    xline(tau_bi, '--'); yline(y_bi, '--');

    xlabel('\tau (s)');
    ylabel('Normal Allan Deviation (m/s^2)');
    title('Accelerometer');
    legend({'X','Y','Z','Average X-Y-Z','Bias instability','White noise'}, 'Location','southeast');

    g0 = 9.80665;
    sigma_rw_mgal_sqrth = sigma_rw * 60;

    txt = sprintf([ ...
        'sigma_bi = %.4g m/s^2 = %.2f mGal\n' ...
        'sigma_rw = %.4g m/s/sqrt(h)' ], ...
        sigma_bi, sigma_bi/1e-5, sigma_rw_mgal_sqrth);

    text(0.03, 0.08, txt, ...
        'Units','normalized', ...
        'Interpreter','none');


end

function plot_gyr(tau, gx, gy, gz, gavg, sigma_rw_rad, sigma_bi_rad, tau_rw, tau_bi)
    figure('Name','Gyroscope Allan Deviation'); clf;
    loglog(tau, gx, 'LineWidth', 1.0); hold on; grid on;
    loglog(tau, gy, 'LineWidth', 1.0);
    loglog(tau, gz, 'LineWidth', 1.0);
    loglog(tau, gavg,'LineWidth', 1.6);

    xlim([1e-2 1e4]);
    set(gca,'XTick',[1e-2 1e-1 1 10 100 1e3 1e4]);

    ref_rw = sigma_rw_rad ./ sqrt(tau);
    ref_bi = sigma_bi_rad * ones(size(tau));
    loglog(tau, ref_bi, '--', 'LineWidth', 1.2);
    loglog(tau, ref_rw, '--', 'LineWidth', 1.2);

    y_rw = sigma_rw_rad / sqrt(tau_rw);
    y_bi = sigma_bi_rad;
    xline(tau_rw, '--'); yline(y_rw, '--');
    xline(tau_bi, '--'); yline(y_bi, '--');

    xlabel('\tau (s)');
    ylabel('Normal Allan Deviation (rad/s)');
    title('Gyroscope');
    legend({'X','Y','Z','Average X-Y-Z','Bias instability','White noise'}, 'Location','southeast');

    sigma_rw_deg_sqrth = (sigma_rw_rad * 180/pi) * 60;
    sigma_bi_degph    = (sigma_bi_rad * 180/pi) * 3600;

    txt = sprintf([ ...
        'sigma_bi = %.2f deg/h\n' ...
        'sigma_rw = %.2f deg/sqrt(h)' ], ...
        sigma_bi_degph, sigma_rw_deg_sqrth);

    text(0.03, 0.08, txt, ...
        'Units','normalized', ...
        'Interpreter','none');

end

% ===================== estimation helpers =====================
function idx = pick_slope_region(tau, adev, targetSlope, tol, whichSide)
    lt = log10(tau(:));
    la = log10(adev(:));
    s  = diff(la) ./ diff(lt);
    s  = [s(1); s(:)];

    cand = find(abs(s - targetSlope) < tol);
    if isempty(cand), idx = 1:min(10,numel(tau)); return; end

    if strcmpi(whichSide,'smallTau')
        cand = cand(tau(cand) <= prctile(tau, 40));
    end
    if isempty(cand), idx = 1:min(10,numel(tau)); return; end

    blocks = split_into_blocks(cand);
    [~, ib] = max(cellfun(@numel, blocks));
    idx = blocks{ib};

    if numel(idx) < 5
        idx = cand(1:min(10,numel(cand)));
    end
end

function blocks = split_into_blocks(idxs)
    d = diff(idxs);
    cut = [0; find(d>1); numel(idxs)];
    blocks = cell(numel(cut)-1,1);
    for i=1:numel(blocks)
        blocks{i} = idxs(cut(i)+1:cut(i+1));
    end
end

function [sigma_bi, idxMin] = bias_from_min(adev)
    [sigmaMin, idxMin] = min(adev);
    sigma_bi = sigmaMin / 0.664;
end

% ===================== Allan deviation core =====================
function [tau, adev] = allan_adev_overlapping(x, dt)
    x = x(:);
    x = x(~isnan(x));
    N = length(x);
    if N < 500
        error('Data too short for stable Allan analysis. Try longer static data.');
    end

    mMax = floor(N/10);
    mMin = 1;
    nTau = 70;
    mList = unique(max(mMin, round(logspace(log10(mMin), log10(mMax), nTau))));
    cs = [0; cumsum(x)];

    tau  = zeros(size(mList));
    adev = zeros(size(mList));

    for i=1:length(mList)
        m = mList(i);
        if (N - 2*m) < 1
            tau(i) = NaN; adev(i) = NaN; continue;
        end
        y = (cs(1+m:end) - cs(1:end-m)) / m;
        d = y(1+m:end) - y(1:end-m);
        avar = 0.5 * mean(d.^2);
        tau(i)  = m * dt;
        adev(i) = sqrt(avar);
    end

    v = isfinite(tau) & isfinite(adev) & tau>0 & adev>0;
    tau = tau(v);
    adev = adev(v);
end