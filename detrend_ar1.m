% smooth time series at specified interval, subtract from time series, 
% bin accordingly, and use for AR-1 model
% 11/20/17 - ECO

clear all 
close all

%%%%% Input parameters %%%%%

% where to find comma seperated data
ptf = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

% define the hour spacing to use for smoothing, binning, and lags
hrs = 24*28; % number of hours for smoothing and filtering
bin_hrs = 24; % number of hours for binning 
std_int = 24; % define the time interval of interest for the lags

% define which data series to use. Options are 'oith', 'para', or 'egg'
ob_str = 'egg'; % string indicating series that will be 'observed'
pd_str = 'oith'; % string indicating series that will be 'predicted'

%%
%%%%% Read data %%%%%
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
pds = pl_holder.(pd_str);

%%%%% Date and time labels for plotting %%%%%
startDate = datenum('03-11-2015 12:00:00 PM');
endDate = datenum('08-01-2015 12:00:00 AM');
xx = linspace(startDate, endDate, max(size(data)));
str = datestr(xx, 'mm-dd');

%%
%%%%% Make moving average filter for smoothing %%%%%
coeffHrs = ones(1, hrs)/hrs;

%%%%% Filter over defined period %%%%%%
avg_obs_hrs = filter(coeffHrs, 1, obs);
avg_pds_hrs = filter(coeffHrs, 1, pds);

%%%%% Smooth, scale, and bin into desired time period
delta_ods = obs - avg_obs_hrs;
delta_pds = pds - avg_pds_hrs;

% scale before binning
%{
%do_scale = (delta_ods - mean(delta_ods))/std(delta_ods);
%dp_scale = (delta_pds - mean(delta_pds))/std(delta_pds);

%do_bin = movmean(do_scale, bin_hrs, 'omitnan');
%dp_bin = movmean(dp_scale, bin_hrs, 'omitnan');
%}

% scale after binning
do_bin = movmean(delta_ods, bin_hrs, 'omitnan');
dp_bin = movmean(delta_pds, bin_hrs, 'omitnan');

do_bin = (do_bin - mean(do_bin))/std(do_bin);
dp_bin = (dp_bin - mean(dp_bin))/std(dp_bin);

%%%%% Get original data in same format to compare to smoothed %%%%%
o_scale = (obs - mean(obs))/std(obs);
p_scale = (pds - mean(pds))/std(pds);

o_bin = movmean(o_scale, bin_hrs, 'omitnan');
p_bin = movmean(p_scale, bin_hrs, 'omitnan');

% plot scaled time series
figure;
plot(xx, p_bin, xx, dp_bin, 'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', xx)
axis tight
datetick('x', 'mm-dd', 'keeplimits')
grid on
title(['Filtered at ', num2str(hrs), ' hour interval'], 'FontSize', 14)
legend([pd_str, ' scaled and binned'], [pd_str, ' smoothed, scaled, and binned'])
xlabel('Time', 'FontSize', 12)
ylabel('scaled number of ROIs', 'FontSize', 12)

figure;
plot(xx, o_bin, xx, do_bin, 'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', xx)
axis tight
datetick('x', 'mm-dd', 'keeplimits')
grid on
title(['Filtered at ', num2str(hrs), ' hour interval'], 'FontSize', 14)
legend([ob_str, ' scaled and binned'], [ob_str, ' smoothed, scaled, and binned'])
xlabel('Time', 'FontSize', 12)
ylabel('scaled number of ROIs', 'FontSize', 12)

%%
%%%%% Run AR-1 model %%%%%

% define lags and pre-allocate space for r2 and significance metrics  
dt = bin_hrs/std_int;
tau = 1:(dt/4):(5/dt);
r_ = zeros(size(tau));
p_ = zeros(size(tau));

% do the prediction
for jj = 1:max(size(tau))
    preds = simple_dynamic_link(do_bin, dt, tau(jj));
    [r_temp, p_temp] = corrcoef(preds, dp_bin); % get the r^2 value
    r_(jj) = r_temp(1,2);
    p_(jj) = p_temp(1,2);
    if jj > 1
        if r_(jj) > r_(jj-1)
            pred_best = preds;
        end
    else
        pred_best = preds;
    end
    clear temp
end

% predicted and actual
figure
titleout2 = [pd_str, ' driven by ', ob_str, ' with ', num2str(hrs), ...
    ' hr detrend: actual data and best AR-1 prediction'];
plot(xx, dp_bin, xx, pred_best, 'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', xx)
axis tight
datetick('x', 'mm-dd', 'keeplimits')
grid on
title(titleout2, 'FontSize', 14)
legend('Measured', 'Predicted')
xlabel('Time', 'FontSize', 12)
ylabel('Scaled Number of ROIs', 'FontSize', 12)

% r^2 series
figure
titleout = [pd_str, ' driven by ', ob_str, ' with ', num2str(hrs), ...
    ' hr detrend; r^2 values'];
plot(tau, r_, '-o', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title(titleout, 'FontSize', 14)
xlabel('time lag', 'FontSize', 12)
ylabel('r^2', 'FontSize', 12)

%%%%% Print the maximum r and associated tau %%%%%
[xx, ii] = max(r_);

sprintf('r^2 max = %0.5g \n tau = %0.5g', xx, tau(ii)) 