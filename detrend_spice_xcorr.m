clear all

%%%%% Input parameters %%%%%

% where to find comma seperated data
ptf = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

% define the hour spacing to use for smoothing, binning, and lags
hrs = 24*7; % number of hours for smoothing and filtering
bin_hrs = 24; % number of hours for binning 
std_int = 24; % define the time interval of interest for the lags

% define which data series to use. Options are 'oith', 'para', or 'egg'
ob_str = 'oith'; % string indicating series that will be 'observed'

%%
%%%%% Read oithona data %%%%%
data = dlmread(ptf, ',');

% get first and last nonzero indicies
nz_front = find(data(:, 2), 1, 'first');
nz_end = find(data(:, 2), 1, 'last');

% pull out the relevent data
data = data(nz_front:nz_end, :); % hard coded to truncate zeros
pl_holder.oith = data(1:end, 2) - data(1:end, 3) + data(1:end, 4);
pl_holder.para = data(1:end, 8) - data(1:end, 9) + data(1:end, 10);
pl_holder.egg = data(1:end, 5) - data(1:end, 6) + data(1:end, 7);

% get the appropriate data
obs = pl_holder.(ob_str);

%% 
%%%%% read data and smooth %%%%%
ff = readtable('/Users/Orenstein/Documents/documents/SPC/oithona_project/autoss_8172_d7f5_b983.csv');

ff = ff(1:34345, :); % just select till august

% smooth over 24 hour periods. Data collected at 6 minute intervals
% multiply bin_hrs by 10 since sampling interval is every 6 minutes
ff.salinity_1e_3_ = movmean(ff.salinity_1e_3_, bin_hrs*10, 'omitnan');
ff.temperature_celsius_ = movmean(ff.temperature_celsius_, bin_hrs*10, 'omitnan');
ff.pressure_dbar_ = movmean(ff.pressure_dbar_, bin_hrs*10, 'omitnan');

% get rid of any bad data points according to data QC flag
p_flag = ff.pressure_flagPrimary==1;
s_flag = ff.salinity_flagPrimary==1;
t_flag = ff.temperature_flagPrimary==1;

flags = p_flag + s_flag + t_flag;
flags = (flags == 3);

xx = ff(flags, :);

% get rid of any questionable salinity points outside of 2 STDs
sd_s = std(xx.salinity_1e_3_);
mu_s = mean(xx.salinity_1e_3_);

under = xx.salinity_1e_3_ > (mu_s - 2*sd_s);
over = xx.salinity_1e_3_ < (mu_s + 2*sd_s);

ou = under + over;
ou = (ou == 2);

xx = xx(ou, :);

clear -VARS under over ou

% get rid of any questionable temperature points outside of 2 STDs
sd_t = std(xx.temperature_celsius_);
mu_t = mean(xx.temperature_celsius_);

under = xx.temperature_celsius_ > (mu_t - 2*sd_t);
over = xx.temperature_celsius_ < (mu_t + 2*sd_t);

ou = under + over;
ou = (ou == 2);

xx = xx(ou, :);

clear -VARS under over ou

% get rid of any questionable pressure points outside of 2 STDs
sd_p = std(xx.pressure_dbar_);
mu_p = mean(xx.pressure_dbar_);

under = xx.pressure_dbar_ > (mu_p - 2*sd_p);
over = xx.pressure_dbar_ < (mu_p + 2*sd_p);

ou = under + over;
ou = (ou == 2);

xx = xx(ou, :);

clear -VARS under over ou

times = xx.time_UTC_;
out = xx{:, {'pressure_dbar_', 'temperature_celsius_', 'salinity_1e_3_'}};

% get date stuff for plotting
startDate = datenum('03-11-2015 12:00:00 PM');
endDate = datenum('08-01-2015 12:00:00 AM');
dd = linspace(startDate, endDate, max(size(xx)));
str = datestr(dd, 'mm-dd');
dd2 = linspace(startDate, endDate, max(size(obs)));

%% GSW spice via TEOS-10
s_r = gsw_SR_from_SP(out(:,3));
c_t = gsw_CT_from_t(s_r, out(:,2), out(:, 1));

gsw_spice = gsw_spiciness0(s_r, c_t);

%%
%%%%% Make moving average filter for smoothing %%%%%
coeffHrs = ones(1, hrs)/hrs;

%%%%% Filter over defined period %%%%%%
avg_obs_hrs = filter(coeffHrs, 1, obs);
avg_spice_hrs = filter(coeffHrs, 1, gsw_spice);

%%%%% Smooth, scale, and bin into desired time period
delta_obs = obs - avg_obs_hrs;
delta_spice = gsw_spice - avg_spice_hrs;

do_scale = (delta_obs - mean(delta_obs))/std(delta_obs);
dsp_scale = (delta_spice - mean(delta_spice))/std(delta_spice);

% try with out scaling (uncomment to use)
%dsp_scale = delta_spice;
%do_scale = movmean(delta_obs, bin_hrs, 'omitnan');

do_bin = movmean(do_scale, bin_hrs, 'omitnan');
%dp_bin = movmean(dp_scale, bin_hrs*10, 'omitnan');

%%
%%%%%% Compute xcorr %%%%%%
% interpolate do_bin to same length as spice
o_interp = interp1(dd2, do_bin, dd); % linear interpolation

[C, lags] = xcorr(o_interp, dsp_scale, 'coeff');

figure;
plot(lags, C, 'LineWidth', 2)
grid on
xlabel('Lags', 'FontSize', 12)
ylabel('Correclation Coefficent', 'FontSize', 12)
set(gca, 'FontSize', 12)

[mx, ii] = max(C);
t_str = {['xcorr with spice lagging ', ob_str, ' detrend at ', num2str(hrs), ' hours'], ...
    sprintf('r^2 max = %0.5g at tau = %0.5g', mx, lags(ii))}; 
title(t_str, 'FontSize', 14)
