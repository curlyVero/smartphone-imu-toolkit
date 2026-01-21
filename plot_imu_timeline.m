% =========================================================================
% plot_imu_timeline - Plot IMU sampling instants for ACC and GYR streams.
%
% DESCRIPTION:
%   Given accelerometer and gyroscope timestamps (UTC milliseconds), this
%   function visualizes sampling instants to inspect alignment, jitter, and
%   missing samples. Gyro timestamps are also rendered as vertical dashed
%   reference lines.
%
% USAGE:
%   plot_imu_timeline(accUtcMs, gyrUtcMs, rate_acc, rate_gyr);
%   plot_imu_timeline(accUtcMs, gyrUtcMs, rate_acc, rate_gyr, true);   % middle 10 s (default)
%   plot_imu_timeline(accUtcMs, gyrUtcMs, rate_acc, rate_gyr, false);  % full duration
%
% INPUTS:
%   accUtcMs, gyrUtcMs : int64/float arrays of UTC milliseconds.
%   rate_acc, rate_gyr : structs returned by imu_rate_stats (fields: f, dt_*).
%   showMid10s         : logical, whether to display only the middle 10 s.
%
% OUTPUT:
%   A figure.
%
% =========================================================================

function plot_imu_timeline(accUtcMs, gyrUtcMs, rate_acc, rate_gyr, showMid10s)
% Plot IMU sampling instants (ACC vs GYR).

    if nargin < 5 || isempty(showMid10s)
        showMid10s = true;
    end

    % ---- time alignment ----
    t0 = min([accUtcMs(:); gyrUtcMs(:)]);
    t_acc = double(accUtcMs - t0) / 1000.0;
    t_gyr = double(gyrUtcMs - t0) / 1000.0;

    % ---- optional: keep only middle 10 seconds ----
    if showMid10s
        t_all = sort([t_acc(:); t_gyr(:)]);
        t_mid = t_all(round(numel(t_all)/2));
        t_start = t_mid - 5.0;
        t_end   = t_mid + 5.0;

        idx_acc = (t_acc >= t_start) & (t_acc <= t_end);
        idx_gyr = (t_gyr >= t_start) & (t_gyr <= t_end);

        t_acc = t_acc(idx_acc);
        t_gyr = t_gyr(idx_gyr);
    end

    % ---- vertical positions ----
    y_acc = 0.45;
    y_gyr = 0.55;

    figure; hold on;
    grid off;
    box on;

    % ---- axis range ----
    xMin = min([t_acc(:); t_gyr(:)]);
    xMax = max([t_acc(:); t_gyr(:)]);

    % ---- gyro vertical dashed lines ----
    xLines = [t_gyr.'; t_gyr.'];
    yLines = [zeros(1, numel(t_gyr)); ones(1, numel(t_gyr))];
    h = line(xLines, yLines);
    set(h, 'LineStyle', '--', 'LineWidth', 0.6, 'Color', [0.6 0.6 0.6]);

    % ---- accelerometer samples ----
    plot(t_acc, y_acc * ones(size(t_acc)), '^', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', [0.2 0.4 0.8], ...
        'MarkerEdgeColor', [0.2 0.4 0.8]);

    % ---- gyroscope samples ----
    plot(t_gyr, y_gyr * ones(size(t_gyr)), 'o', ...
        'MarkerSize', 5, ...
        'MarkerFaceColor', [0.85 0.33 0.1], ...
        'MarkerEdgeColor', [0.85 0.33 0.1]);

    % ---- axes formatting ----
    xlim([xMin xMax]);
    ylim([0 1]);
    yticks([y_acc y_gyr]);
    yticklabels({'ACC', 'GYR'});

    xlabel('Time since start (s)');
    if showMid10s
        title('IMU Sampling Timeline (Middle 10 s)');
    else
        title('IMU Sampling Timeline (Full Duration)');
    end

    legend({ ...
        'GYR sample instants', ...
        sprintf('ACC (%.2f Hz)', rate_acc.f), ...
        sprintf('GYR (%.2f Hz)', rate_gyr.f)}, ...
        'Location', 'northoutside', ...
        'Orientation', 'horizontal');
end
