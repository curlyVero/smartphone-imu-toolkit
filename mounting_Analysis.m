% =========================================================================
% mounting_Analysis - Estimate smartphone IMU mounting angles in a vehicle.
%
% DESCRIPTION:
%   This tool estimates the mounting (installation) angles between the
%   smartphone body frame and the vehicle frame in two phases:
%
%   Phase 1 (static leveling):
%     Estimate roll/pitch by aligning the mean accelerometer vector during
%     detected static segments to gravity. The optimization minimizes the
%     horizontal (X/Y) components of the rotated gravity vector.
%
%   Phase 2 (yaw heading):
%     After roll/pitch leveling, estimate yaw by maximizing the correlation
%     between lateral acceleration and yaw rate during dynamic turning
%     segments. A simple grid search is used.
%
% USAGE:
%   result = mounting_Analysis('imu_preprocess/example.20IMU');
%   result = mounting_Analysis('imu_preprocess/example.20IMU', 'Plot', true);
%
% INPUT:
%   imuFile: Preprocessed IMU file (CSV/text) with columns:
%     [t_ms, AccX, AccY, AccZ, GyroX, GyroY, GyroZ]
%   Gyro is assumed rad/s. Acc is m/s^2.
%
% NAME-VALUE OPTIONS:
%   'Gravity'          : gravity magnitude in m/s^2 (default 9.80665)
%   'Plot'             : whether to generate diagnostic figures (default true)
%   'WinSize'          : window length for static detection (default 50)
%   'StaticStdThresh'  : std threshold on |acc| for static detection (default 0.25)
%   'AccMoveThresh'    : dynamic threshold on |acc_dyn| (default 0.2)
%   'YawRateThresh'    : dynamic threshold on |yaw_rate| (default 0.02)
%   'YawGridStepDeg'   : yaw grid step in degrees (default 0.1)
%
% OUTPUT:
%   result: struct with fields (degrees unless noted):
%     roll_deg, pitch_deg, yaw_deg, max_corr
%     static_mean_b (1x3) mean acc in body frame (m/s^2)
%     n_static, n_total, dt, fs
%
% NOTES:
%   - If Aerospace Toolbox is unavailable, a lightweight XYZ Euler-to-DCM
%     implementation is used internally.
%
% =========================================================================
function result = mounting_Analysis(imuFile, varargin)

    if nargin < 1 || isempty(imuFile)
        error('imuFile is required. Example: result = mounting_Analysis(''imu_preprocess/example.20IMU'');');
    end

    p = inputParser;
    p.addRequired('imuFile', @(x) ischar(x) || isstring(x));
    p.addParameter('Gravity', 9.80665, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('Plot', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('WinSize', 50, @(x) isnumeric(x) && isscalar(x) && x >= 5);
    p.addParameter('StaticStdThresh', 0.25, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('AccMoveThresh', 0.2, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('YawRateThresh', 0.02, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('YawGridStepDeg', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.parse(imuFile, varargin{:});

    imuFile = char(p.Results.imuFile);
    GRAVITY_MAG = p.Results.Gravity;
    doPlot = logical(p.Results.Plot);
    win_size = p.Results.WinSize;
    static_std_th = p.Results.StaticStdThresh;
    acc_move_th = p.Results.AccMoveThresh;
    yaw_rate_th = p.Results.YawRateThresh;
    yaw_step = p.Results.YawGridStepDeg;

    fprintf('Loading file: %s ...\n', imuFile);

    % ---- load data ----
    try
        opts = detectImportOptions(imuFile, 'FileType', 'text');
        opts.DataLines = [2 Inf];
        opts.VariableNamingRule = 'preserve';
        raw_table = readtable(imuFile, opts);

        if ismember('AccX', raw_table.Properties.VariableNames)
            acc_b  = [raw_table.AccX,  raw_table.AccY,  raw_table.AccZ];
            gyro_b = [raw_table.GyroX, raw_table.GyroY, raw_table.GyroZ];
            t_ms   = raw_table{:,1};
        else
            data_mat = table2array(raw_table);
            t_ms   = data_mat(:,1);
            acc_b  = data_mat(:,2:4);
            gyro_b = data_mat(:,5:7);
        end
    catch ME
        error('File load error: %s', ME.message);
    end

    t = (t_ms - t_ms(1)) * 1e-3;
    dt = median(diff(t));
    fs = 1/dt;

    fprintf('Data loaded. Points: %d (fs ~ %.3f Hz)\n', size(acc_b,1), fs);

    % ---- phase 1: roll/pitch from static segments ----
    acc_norm = vecnorm(acc_b, 2, 2);
    acc_std  = movstd(acc_norm, win_size);
    static_mask = acc_std < static_std_th;

    if sum(static_mask) < 10
        warning('Not enough static samples were detected. Using mean of all data.');
        static_mean = mean(acc_b, 1);
    else
        static_mean = mean(acc_b(static_mask, :), 1);
    end

    fprintf('Static mean acc (body frame): [%.4f, %.4f, %.4f]\n', static_mean);
    n_static = sum(static_mask);

    % Optimize roll/pitch so that rotated gravity has minimal horizontal component
    target_xy = [0; 0];
    S = [1, 0, 0; 0, 1, 0];

    cost_func_rp = @(x) norm( ...
        S * (angle2dcm_xyz(x(1), -x(2), 0) * static_mean') - target_xy ...
        )^2;

    phi_init   = atan2(static_mean(2), -static_mean(3));
    theta_init = atan2(static_mean(1), sqrt(static_mean(2)^2 + static_mean(3)^2));
    x0 = [phi_init, theta_init];

    options = optimset('Display', 'off', 'TolX', 1e-8);
    x_opt = fminsearch(@(x) cost_func_rp([x(1), x(2), 0]), x0, options);

    best_roll  = x_opt(1);
    best_pitch = x_opt(2);

    fprintf('--------------------------------------\n');
    fprintf('Phase 1: Roll/Pitch optimized\n');
    fprintf('  Roll  : %.4f deg\n', rad2deg(best_roll));
    fprintf('  Pitch : %.4f deg\n', rad2deg(best_pitch));

    R_level = angle2dcm_xyz(best_roll, -best_pitch, 0);

    % ---- phase 2: yaw from correlation between lateral acc and yaw rate ----
    acc_level  = (R_level * acc_b')';
    gyro_level = (R_level * gyro_b')';

    % Remove gravity in the leveled frame (Down is -Z in an FRD convention)
    acc_dyn_level = acc_level - [0, 0, -GRAVITY_MAG];

    yaw_rate_meas = gyro_level(:, 3);
    acc_dyn_norm  = vecnorm(acc_dyn_level, 2, 2);

    corr_mask = (acc_dyn_norm > acc_move_th) & (abs(yaw_rate_meas) > yaw_rate_th);
    valid_dyn  = acc_dyn_level(corr_mask, :);
    valid_rate = yaw_rate_meas(corr_mask);

    if numel(valid_rate) < 20
        warning('Not enough dynamic samples for yaw estimation. Yaw may be unreliable.');
    end

    yaw_grid = -180:yaw_step:180;
    corrs = zeros(size(yaw_grid));

    for i = 1:length(yaw_grid)
        psi = deg2rad(yaw_grid(i));
        % Rotate leveled acceleration around Z (heading)
        lat_acc_test = -valid_dyn(:, 1) * sin(psi) + valid_dyn(:, 2) * cos(psi);

        if std(lat_acc_test) > 1e-6 && std(valid_rate) > 1e-6
            c = corrcoef(lat_acc_test, valid_rate);
            corrs(i) = c(1, 2);
        else
            corrs(i) = 0;
        end
    end

    [max_corr, idx] = max(corrs);
    best_yaw = yaw_grid(idx);

    fprintf('--------------------------------------\n');
    fprintf('Phase 2: Yaw optimized (target: high positive correlation)\n');
    fprintf('  Yaw     : %.2f deg\n', best_yaw);
    fprintf('  Max Corr: %.4f\n', max_corr);

    % ---- pack results ----
    result = struct();
    result.roll_deg  = rad2deg(best_roll);
    result.pitch_deg = rad2deg(best_pitch);
    result.yaw_deg   = best_yaw;
    result.max_corr  = max_corr;
    result.static_mean_b = static_mean;
    result.n_static = n_static;
    result.n_total  = size(acc_b,1);
    result.dt = dt;
    result.fs = fs;

    % ---- diagnostics plot ----
    if doPlot
        R_final = angle2dcm_xyz(best_roll, -best_pitch, deg2rad(best_yaw));
        acc_frd = (R_final * acc_b')';
        acc_frd_dyn = acc_frd - [0, 0, -GRAVITY_MAG];

        figure('Name', 'Mounting Analysis', 'Color', 'w', 'Position', [100, 100, 1400, 900]);

        % Time-series sign check
        subplot(2, 2, 1);
        plot_range = 1:length(yaw_rate_meas);
        if length(plot_range) > 2000
            [~, center] = max(abs(yaw_rate_meas));
            plot_range = max(1, center-1000) : min(length(yaw_rate_meas), center+1000);
        end
        yyaxis left;
        plot(plot_range, acc_frd_dyn(plot_range, 2), 'LineWidth', 1.5);
        ylabel('Lateral Acc (m/s^2)');
        yyaxis right;
        plot(plot_range, yaw_rate_meas(plot_range), '--', 'LineWidth', 1.5);
        ylabel('Yaw Rate (rad/s)');
        title(sprintf('Sign Check: LatAcc vs YawRate (corr = %.2f)', max_corr));
        legend('Lat Acc', 'Yaw Rate');
        grid on;

        % Scatter in F-R plane
        subplot(2, 2, 2);
        high_dyn_mask = abs(yaw_rate_meas) > yaw_rate_th;
        if sum(high_dyn_mask) > 10
            scatter(acc_frd_dyn(high_dyn_mask, 1), acc_frd_dyn(high_dyn_mask, 2), ...
                20, abs(yaw_rate_meas(high_dyn_mask)), 'filled');
            colorbar;
            title('Forward vs Lateral Acc (colored by |YawRate|)');
            xlabel('Forward Acc (m/s^2)');
            ylabel('Lateral Acc (m/s^2)');
            grid on; axis equal;
            xline(0, '--'); yline(0, '--');

            pfit = polyfit(acc_frd_dyn(high_dyn_mask, 1), acc_frd_dyn(high_dyn_mask, 2), 1);
            ref_x = linspace(min(acc_frd_dyn(high_dyn_mask, 1)), max(acc_frd_dyn(high_dyn_mask, 1)), 100);
            hold on; plot(ref_x, polyval(pfit, ref_x), '-', 'LineWidth', 2);
            text(min(xlim)+0.5, max(ylim)-0.5, sprintf('Slope: %.3f', pfit(1)), 'BackgroundColor', 'w');
        else
            text(0.5, 0.5, 'Not enough dynamic data', 'HorizontalAlignment', 'center');
        end

        % Correlation vs yaw grid
        subplot(2, 2, 3);
        plot(yaw_grid, corrs, 'LineWidth', 1.5);
        hold on; xline(best_yaw, '--');
        title('Yaw Grid Search');
        ylabel('Correlation coefficient'); xlabel('Yaw angle (deg)');
        grid on;
        text(best_yaw, max_corr, sprintf('  best: %.2f deg', best_yaw));

        % Straight-motion scatter (yaw_rate ~ 0)
        subplot(2, 2, 4);
        straight_mask = abs(yaw_rate_meas) < max(0.05, yaw_rate_th) & vecnorm(acc_frd_dyn, 2, 2) > 0.5;
        if sum(straight_mask) > 10
            scatter(acc_frd_dyn(straight_mask, 1), acc_frd_dyn(straight_mask, 2), ...
                15, 'filled', 'MarkerFaceAlpha', 0.5);
            title('Forward vs Lateral Acc (near-straight motion)');
            xlabel('Forward Acc'); ylabel('Lateral Acc');
            grid on; axis equal;
            yline(0, '-', 'LineWidth', 2);
        else
            text(0.5, 0.5, 'No near-straight segments found', 'HorizontalAlignment', 'center');
        end
    end
end

% ---- local helper: XYZ Euler angles to DCM ----
function R = angle2dcm_xyz(phi, theta, psi)
% angle2dcm_xyz - Minimal replacement for angle2dcm(phi,theta,psi,'XYZ').
%
% Rotation order: X then Y then Z.
    cx = cos(phi);  sx = sin(phi);
    cy = cos(theta); sy = sin(theta);
    cz = cos(psi);  sz = sin(psi);

    Rx = [1 0 0; 0 cx -sx; 0 sx cx];
    Ry = [cy 0 sy; 0 1 0; -sy 0 cy];
    Rz = [cz -sz 0; sz cz 0; 0 0 1];

    R = Rx * Ry * Rz;
end
