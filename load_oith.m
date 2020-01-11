function [data, times] = load_oith(path_to_file, start_time, end_time)

%%%%% Load data and print out appriate time scale %%%%%
% 
% Load oithona data assuming csv file with rows as time and cols as 
%
% [timestamp, oith_count, oith_false_positive, oith_mislabeled, ...
%   egg_count, egg_false_positive, egg_mislabeled, ...
%       para_count, para_false_positive, para_mislabeled]
%
%   :param path_to_file: absolute path to csv file [string]
%   :param start_time: start date as YYYY-MM-DD HH:MM:SS [string]
%   :param end_time: end date as YYYY-MM-DD HH:MM:SS [string]
%   :return data: data with corrected counts [cell array]
%   :return times: time information for plotting [cell array]

input = dlmread(path_to_file, ',');

% get first and last nonzero indicies
nz_front = find(input(:, 2), 1, 'first');
nz_end = find(input(:, 2), 1, 'last');

% pull out the relevent data
input = input(nz_front:nz_end, :); % hard coded to truncate zeros
data.oith = input(1:end, 2) - input(1:end, 3) + input(1:end, 4);
data.para = input(1:end, 8) - input(1:end, 9) + input(1:end, 10);
data.egg = input(1:end, 5) - input(1:end, 6) + input(1:end, 7);

%%%%% Date and time labels for plotting %%%%%
startDate = datenum(start_time);
endDate = datenum(end_time);
times.xx = linspace(startDate, endDate, max(size(input)));
times.string = datestr(times.xx, 'mm-dd');

