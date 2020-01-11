% smooth time series at specified interval, subtract from time series, 
% bin accordingly, and use for two tau model

% not working quite right. axes for mesh plot are wacky. 
% r^2 values not substantially better for any combination of observed and
% predicited. 
% 11/30/17 - ECO

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
ob_str = 'para'; % string indicating series that will be 'observed'
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

do_scale = (delta_ods - mean(delta_ods))/std(delta_ods);
dp_scale = (delta_pds - mean(delta_pds))/std(delta_pds);

do_bin = movmean(do_scale, bin_hrs, 'omitnan');
dp_bin = movmean(dp_scale, bin_hrs, 'omitnan');

%%
%%%%% Run AR-1 model %%%%%

% define lags and pre-allocate space for r2 and significance metrics  
dt = bin_hrs/std_int;
[tau1, tau2] = meshgrid((1/dt):(dt/4):(5/dt));

rr = zeros(numel(tau1),1);
% do the prediction
for jj = 1:numel(tau1)
    preds = dynamic_link_posneg(do_bin, dt, tau1(jj), tau2(jj));
    temp = corrcoef(preds, dp_bin); % get the r^2 value
    rr(jj) = temp(1, 2);
    if jj > 1
        if rr(jj) > rr(jj-1)
            pred_best = preds;
        end
    else
        pred_best = preds;
    end
    clearvars temp
end

%reshape rr and save it out
rr_out = reshape(rr, size(tau1, 1), size(tau1, 2));

% predicted and actual
figure
titleout2 = [pd_str, ' driven by ', ob_str, ' with ', num2str(hrs), ...
    ' hr detrend: actual data and best two tau prediction'];
plot(xx, dp_bin, xx, pred_best, 'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', xx)
axis tight
datetick('x', 'mm-dd', 'keeplimits')
set(gca, 'FontSize', 12)
grid on
title(titleout2, 'FontSize', 14)
legend('Measured', 'Predicted')
xlabel('Time', 'FontSize', 12)
ylabel('Scaled Number of ROIs', 'FontSize', 12)

% r^2 series
figure
titleout = [pd_str, ' driven by ', ob_str, ' with ', num2str(hrs), ...
    ' hr detrend; r^2 values'];
surf(tau1, tau2, rr_out)
set(gca, 'FontSize', 12)
title(titleout, 'FontSize', 14)
xlabel('tau1 (pos time lag)', 'FontSize', 12)
ylabel('tau2 (neg time lag)', 'FontSize', 12)
zlabel('r^2', 'FontSize', 12)
grid on


%%%%% Print the maximum r and associated tau %%%%%
[xx, ind] = max(rr_out(:));
[ii, jj] = ind2sub(size(rr_out), ind);

sprintf('r^2 max = %0.5g \n tau1 = %0.5g \n tau2 = %0.5g',...
    xx, tau1(ii), tau2(jj)) 
