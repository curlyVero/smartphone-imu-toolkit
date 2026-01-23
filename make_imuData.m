% =========================================================================
% make_imuData - Merge separate smartphone accelerometer and gyroscope logs
%               into a single "Uncal" CSV (UncalAccel / UncalGyro).
%
% DESCRIPTION:
%   This utility reads two sensor text files (accelerometer and gyroscope)
%   exported separately, converts their timestamps to UTC milliseconds
%   (Unix epoch), and writes a merged CSV in the Android "Uncal" style:
%     MessageType,utcTimeMillis,MeasurementX,MeasurementY,MeasurementZ,BiasX,BiasY,BiasZ
%   The output is sorted by utcTimeMillis and can be used by downstream
%   tools in this repository (e.g., plotting and preprocessing).
%
% USAGE:
%   make_imuData('acc.txt', 'gyr.txt', 'imuData.csv');
%
% INPUT FORMAT (accFile / gyrFile):
%   Whitespace-separated numeric columns (header lines are allowed):
%     year month day hour minute second x y z [biasx biasy biasz]
%   - If bias columns are missing, biases are filled with 0.
%   - Units:
%       Accelerometer: m/s^2
%       Gyroscope:     rad/s
%
% OUTPUT:
%   outFile: CSV with header
%     MessageType,utcTimeMillis,MeasurementX,MeasurementY,MeasurementZ,BiasX,BiasY,BiasZ
%
% DEPENDENCIES:
%   plot_imu_timeline.m (optional visualization)
%
% =========================================================================

function make_imuData(accFile, gyrFile, outFile)
% Input (whitespace separated, header optional):
%  year month day hour minute second x y z [biasx biasy biasz]
%
% Output CSV with header:
%  MessageType,utcTimeMillis,MeasurementX,MeasurementY,MeasurementZ,BiasX,BiasY,BiasZ

    if nargin < 3
        error('Usage: make_imuData(accFile, gyrFile, outFile)');
    end

    % ===== Read =====
    acc = read_numeric_sensor(accFile);
    gyr = read_numeric_sensor(gyrFile);

    % ===== Time to UTC millis =====
    accUtcMs = to_utc_millis(acc);
    gyrUtcMs = to_utc_millis(gyr);

    % ===== Sampling rate stats (computed once) =====
    rate_acc = imu_rate_stats(accUtcMs);
    rate_gyr = imu_rate_stats(gyrUtcMs);

    print_rate("ACC", rate_acc);
    print_rate("GYR", rate_gyr);

    % ===== Plot timeline =====
    % plot_imu_timeline(accUtcMs, gyrUtcMs, rate_acc, rate_gyr);

    % ===== Pack output rows =====
    accRows = pack_rows("UncalAccel", acc, accUtcMs);
    gyrRows = pack_rows("UncalGyro",  gyr, gyrUtcMs);

    allRows = [accRows; gyrRows];
    allRows = sortrows(allRows, "utcTimeMillis");

     % ===== Write CSV =====
    fid = fopen(outFile, 'w');
    assert(fid > 0, "Cannot open output file: %s", outFile);
    N = height(allRows);
    fprintf('[%s] Writing CSV (%d rows) to %s ...\n', mfilename, N, outFile);
    % Print progress every 2 seconds
    printIntervalSec = 2;
    lastPrint = tic;

    tWrite = tic;

    fprintf(fid, 'MessageType,utcTimeMillis,MeasurementX,MeasurementY,MeasurementZ,BiasX,BiasY,BiasZ\n');
    for i = 1:height(allRows)
        fprintf(fid, '%s,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', ...
            allRows.MessageType(i), ...
            allRows.utcTimeMillis(i), ...
            allRows.MeasurementX(i), allRows.MeasurementY(i), allRows.MeasurementZ(i), ...
            allRows.BiasX(i), allRows.BiasY(i), allRows.BiasZ(i));
        % ----- Progress printing (every 2 seconds) -----
        if toc(lastPrint) >= printIntervalSec
            pct = (double(i) / double(N) * 100);
            fprintf('[%s] Writing CSV ... ~%.2f%% (i=%d/%d)\n', mfilename, pct, i, N);
            lastPrint = tic;
        end
    end
    fclose(fid);

    fprintf('[%s] Writing CSV ... 100%% (done, %.2f s)\n', mfilename, toc(tWrite));
    fprintf('Wrote %d IMU records to %s\n', N, outFile);

end

% =========================
% Readers
% =========================
function T = read_numeric_sensor(fname)
% Robust numeric reader (ignores header lines automatically)
% Accepts 9 columns (no bias) or 12 columns (with bias)
%
% Progress:
%   Prints progress to command line every 2 seconds based on ftell/bytes.
%   Also prints valid row count to confirm it is still running.

    fid = fopen(fname, 'r');
    assert(fid > 0, "Cannot open file: %s", fname);

    info = dir(fname);
    totalBytes = max(1, info.bytes);

    % ---------- Preallocation (key speedup) ----------
    blockSize = 200000;           % adjust: 1e5~5e5
    rows = zeros(blockSize, 12);
    n = 0;

    % ---------- Progress print config ----------
    printIntervalSec = 2;
    lastPrint = tic;

    fprintf('[%s] Reading %s ...\n', mfilename, fname);

    try
        while true
            line = fgetl(fid);
            if ~ischar(line), break; end

            line = strtrim(line);
            if isempty(line), continue; end

            % Fast skip for obvious header lines
            c = line(1);
            if ~((c >= '0' && c <= '9') || c == '-' || c == '.')
                continue;
            end

            nums = sscanf(line, '%f');
            if numel(nums) < 9
                continue;
            end

            % Normalize to 12 columns
            if numel(nums) >= 12
                nums = nums(1:12);
            else
                nums = [nums(1:9); 0; 0; 0];
            end

            % Store (preallocated)
            n = n + 1;
            if n > size(rows,1)
                rows(end+1:end+blockSize, :) = 0; %#ok<AGROW>
            end
            rows(n, :) = nums(:).';

            % ----- Progress printing (every 2 seconds) -----
            if toc(lastPrint) >= printIntervalSec
                curBytes = ftell(fid);
                pct = floor(double(curBytes) / double(totalBytes) * 100);
                pct = min(max(pct, 0), 100);
                fprintf('[%s] Reading %s ... ~%d%% (valid rows=%d)\n', ...
                    mfilename, fname, pct, n);
                lastPrint = tic;
            end
        end
    catch ME
        fclose(fid);
        rethrow(ME);
    end

    fclose(fid);

    fprintf('[%s] Reading %s ... 100%% (done, valid rows=%d)\n', mfilename, fname, n);

    assert(n > 0, "No valid numeric data found in: %s", fname);

    rows = rows(1:n, :);

    T = table;
    T.year   = rows(:,1);
    T.month  = rows(:,2);
    T.day    = rows(:,3);
    T.hour   = rows(:,4);
    T.minute = rows(:,5);
    T.second = rows(:,6);
    T.x      = rows(:,7);
    T.y      = rows(:,8);
    T.z      = rows(:,9);
    T.biasx  = rows(:,10);
    T.biasy  = rows(:,11);
    T.biasz  = rows(:,12);
end


function utcMillis = to_utc_millis(T)
% Convert (Y,M,D,h,m,sec) to UTC milliseconds since 1970-01-01

    secInt  = floor(T.second);
    secFrac = T.second - secInt;

    dt = datetime(T.year, T.month, T.day, T.hour, T.minute, secInt, 'TimeZone', 'UTC');
    utcMillis = int64(posixtime(dt) * 1000 + round(secFrac * 1000));
end

% =========================
% Output packer
% =========================
function outT = pack_rows(msgType, T, utcMillis)
    n = height(T);
    outT = table;
    outT.MessageType   = repmat(string(msgType), n, 1);
    outT.utcTimeMillis = utcMillis;
    outT.MeasurementX  = T.x;
    outT.MeasurementY  = T.y;
    outT.MeasurementZ  = T.z;
    outT.BiasX         = T.biasx;
    outT.BiasY         = T.biasy;
    outT.BiasZ         = T.biasz;
end

% =========================
% Sampling rate (single source of truth)
% =========================
function rate = imu_rate_stats(utcMillis)
% Robust sampling rate statistics based on median dt

    rate = struct('f', NaN, 'dt_min', NaN, 'dt_med', NaN, 'dt_max', NaN, 'n', numel(utcMillis));

    if numel(utcMillis) < 2
        return;
    end

    t = double(utcMillis(:)) / 1000.0; % seconds
    dt = diff(t);
    dt = dt(dt > 0); % remove non-increasing

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
        fprintf('%s: insufficient or invalid timestamps (N=%d)\n', tag, rate.n);
        return;
    end

    fprintf('%s rate ~ %.3f Hz (N=%d, dt min/med/max = %.6f / %.6f / %.6f s)\n', ...
        tag, rate.f, rate.n, rate.dt_min, rate.dt_med, rate.dt_max);
end



