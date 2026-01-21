% =========================================================================
% main_generate_kml - Generate a KML trajectory from a navigation text file.
%
% DESCRIPTION:
%   This utility reads a navigation/trajectory text file (e.g., GNSS/INS
%   solution output) that contains a header line with keywords:
%     "GPSTime", "Latitude", "Longitude"
%   followed by numeric data lines. Latitude/longitude are assumed to be
%   stored as degrees/minutes/seconds (DMS) triplets and are converted to
%   decimal degrees. The trajectory is then exported to a KML file using
%   LineKml.m.
%
% USAGE:
%   main_generate_kml('solution.txt', 'trajectory.kml');
%   main_generate_kml('solution.txt', 'trajectory.kml', 'GNSS_INS_Trajectory');
%
% INPUTS:
%   inputFile  : Path to the navigation text file.
%   outputKml  : Output KML path.
%   trackName  : (Optional) Name shown in Google Earth (default: 'Trajectory').
%
% DEPENDENCIES:
%   LineKml.m
%
% =========================================================================
function main_generate_kml(inputFile, outputKml, trackName)

    % Backward-friendly defaults:
    % If called with no arguments, try ./output/nav_antenna.txt under the current folder.
    if nargin < 1 || isempty(inputFile)
        inputFile = fullfile(pwd, 'output', 'nav_antenna.txt');
        fprintf('No inputFile specified. Using default: %s\n', inputFile);
    end
    if nargin < 2 || isempty(outputKml)
        [p, n] = fileparts(inputFile);
        outputKml = fullfile(p, [n, '_trajectory.kml']);
    end
    if nargin < 3 || isempty(trackName)
        trackName = 'Trajectory';
    end

    fprintf('Reading input file: %s\n', inputFile);
    [latitude, longitude, height] = readNavigationData(inputFile);

    if isempty(latitude)
        error('No valid data points were read from: %s', inputFile);
    end

    fprintf('Read %d data points.\n', length(latitude));

    fprintf('Generating KML: %s\n', outputKml);
    fileID = fopen(outputKml, 'w');
    if fileID == -1
        error('Cannot create KML file: %s', outputKml);
    end

    LineKml(fileID, trackName, latitude, longitude, height);

    fclose(fileID);
    fprintf('KML file generated: %s\n', outputKml);
end

function [latitude, longitude, height] = readNavigationData(filename)
% Read a navigation text file and extract latitude/longitude/height.
%
% The function searches for the header line containing "GPSTime", "Latitude"
% and "Longitude", then skips one unit line, and parses subsequent numeric
% lines. Expected columns (at least):
%   [GPSTime, lat_deg, lat_min, lat_sec, lon_deg, lon_min, lon_sec, height, ...]
%
% OUTPUT:
%   latitude/longitude in decimal degrees, height in meters.

    latitude  = [];
    longitude = [];
    height    = [];

    fid = -1;
    try
        fid = fopen(filename, 'r');
        if fid == -1
            error('Cannot open file: %s', filename);
        end

        % Locate the data section
        dataFound = false;
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line) && contains(line, 'GPSTime') && contains(line, 'Latitude') && contains(line, 'Longitude')
                fgetl(fid); % skip the unit line
                dataFound = true;
                break;
            end
        end

        if ~dataFound
            error('Cannot find the beginning of the data section (header not found).');
        end

        tempLat = [];
        tempLon = [];
        tempH   = [];

        while ~feof(fid)
            line = fgetl(fid);
            if ~ischar(line) || isempty(strtrim(line)) || ~contains(line, '.')
                continue;
            end

            data = sscanf(line, ...
                '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', ...
                [1, inf]);

            if numel(data) >= 8
                lat_deg = data(2); lat_min = data(3); lat_sec = data(4);
                lon_deg = data(5); lon_min = data(6); lon_sec = data(7);

                tempLat(end+1) = lat_deg + lat_min/60 + lat_sec/3600; %#ok<AGROW>
                tempLon(end+1) = lon_deg + lon_min/60 + lon_sec/3600; %#ok<AGROW>
                tempH(end+1)   = data(8); %#ok<AGROW>
            end
        end

        fclose(fid);
        fid = -1;

        latitude  = tempLat;
        longitude = tempLon;
        height    = tempH;

    catch ME
        if fid ~= -1
            fclose(fid);
        end
        warning('Error while reading %s: %s', filename, ME.message);
    end
end
