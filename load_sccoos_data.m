function [out] = load_sccoos_data(path_to_file)

%%%%% Load SCCOOS data from csv %%%%%
%
% Assumes can be read in as a table with the standard scoos column names
% Cleans data by removing entries flagged by SCCOOS QC code and removes
% data outside 2 standard deviations.
%
% NOTE: currently have index selecting until the first of august hardcoded
% can make dynamic by editting at line 19 [ECO 12-16-17]
%
%   :param path_to_file: absolute path to csv [string]
%   :return out: cleaned data with cols (pressure, temp, salinity) [array]
%

%% read data and smooth
ff = readtable(path_to_file);

ff = ff(1:34345, :); % just select till august [HARDCODED]

% smooth over 24 hour periods. Data collected at 6 minute intervals
ff.salinity_1e_3_ = movmean(ff.salinity_1e_3_, 240, 'omitnan');
ff.temperature_celsius_ = movmean(ff.temperature_celsius_, 240, 'omitnan');
ff.pressure_dbar_ = movmean(ff.pressure_dbar_, 240, 'omitnan');

%% get rid of any bad data points according to data QC flag
p_flag = ff.pressure_flagPrimary==1;
s_flag = ff.salinity_flagPrimary==1;
t_flag = ff.temperature_flagPrimary==1;

flags = p_flag + s_flag + t_flag;
flags = (flags == 3);

xx = ff(flags, :);

%% get rid of any questionable salinity points outside of 2 STDs
sd_s = std(xx.salinity_1e_3_);
mu_s = mean(xx.salinity_1e_3_);

under = xx.salinity_1e_3_ > (mu_s - 2*sd_s);
over = xx.salinity_1e_3_ < (mu_s + 2*sd_s);

ou = under + over;
ou = (ou == 2);

xx = xx(ou, :);

clear -VARS under over ou

%% get rid of any questionable temperature points outside of 2 STDs
sd_t = std(xx.temperature_celsius_);
mu_t = mean(xx.temperature_celsius_);

under = xx.temperature_celsius_ > (mu_t - 2*sd_t);
over = xx.temperature_celsius_ < (mu_t + 2*sd_t);

ou = under + over;
ou = (ou == 2);

xx = xx(ou, :);

clear -VARS under over ou

%% get rid of any questionable pressure points outside of 2 STDs
sd_p = std(xx.pressure_dbar_);
mu_p = mean(xx.pressure_dbar_);

under = xx.pressure_dbar_ > (mu_p - 2*sd_p);
over = xx.pressure_dbar_ < (mu_p + 2*sd_p);

ou = under + over;
ou = (ou == 2);

xx = xx(ou, :);

clear -VARS under over ou

%times = xx.time_UTC_;
out = xx{:, {'pressure_dbar_', 'temperature_celsius_', 'salinity_1e_3_'}};