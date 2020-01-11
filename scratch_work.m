clear all 
close all

%%%%% Read data %%%%%
data = dlmread('/Users/Orenstein/Desktop/corrected_counts_091817.txt',',');

%%%%% Truncate zeros %%%%%
data = data(26:end-744, :);

%%%%% Make the corrected hourly counts and bin %%%%%
input1 = {'oith', data(1:end-1, 2) - data(1:end-1, 3) + data(1:end-1, 4)};
input2 = {'para', data(1:end-1, 8) - data(1:end-1, 9) + data(1:end-1, 10)};
input3 = {'egg', data(1:end-1, 5) - data(1:end-1, 6) + data(1:end-1, 7)};
time_frames = [1, 6, 24]; % time frames of interest

mm = 'run_mean'; % define method for averaging
%mm = 'block';

dd = bin_pack(time_frames, mm, input1, input2, input3);

%%


oith = data(1:end, 2) - data(1:end, 3) + data(1:end, 4);
oith_scaled = (oith - mean(oith))/std(oith);
oith_fft = fft(oith);

fs = 1;
nn = length(oith);

%fnorm = 

%%
clear all 
close all

%%%%% Read data %%%%%
data = dlmread('/Users/Orenstein/Desktop/corrected_counts_091817.txt',',');

%%%%% Truncate zeros, correct counts, and scale %%%%%
data = data(26:end-744, :);
oith = data(1:end, 2) - data(1:end, 3) + data(1:end, 4);
oith_scaled = (oith - mean(oith))/std(oith);

para = data(1:end, 8) - data(1:end, 9) + data(1:end, 10);
para_scaled = (para - mean(para))/std(para);

%%%%% Date and time labels for plotting %%%%%
startDate = datenum('03-11-2015 12:00:00 PM');
endDate = datenum('08-01-2015 12:00:00 AM');
xx = linspace(startDate, endDate, max(size(data)));
str = datestr(xx, 'mm-dd');

%%%%% Make moving average filter for smoothing %%%%%
hrs = 672;
coeffHrs = ones(1, hrs)/hrs;

%%%%% Filter over defined period %%%%%%
avg_oith_hrs = filter(coeffHrs, 1, oith);
avg_para_hrs = filter(coeffHrs, 1, para);

%%%%% Smooth and bin into desired time period
delta_oith = oith - avg_oith_hrs;
delta_para = para - avg_para_hrs;

do_scale = (delta_oith - mean(delta_oith))/std(delta_oith);
dp_scale = (delta_para - mean(delta_para))/std(delta_para);

do_24 = movmean(do_scale, 24, 'omitnan');
dp_24 = movmean(dp_scale, 24, 'omitnan');

figure
plot(xx, do_24, xx, dp_24)
labs = get(gca,'XLim');
set(gca, 'XTick', xx)
axis tight
datetick('x', 'mm-dd', 'keeplimits')
grid on
title('Filtered at 28 day interval', 'FontSize', 14)
legend('Oithona', 'Paradinium')
xlabel('Time', 'FontSize', 12)
ylabel('scaled number of ROIs', 'FontSize', 12)

dt = 1;
%tau = (1/dt):(dt/4):(5/dt);
tau = 1:(dt/4):(5/dt);
r_ = zeros(size(tau));
p_ = zeros(size(tau));

% do the prediction
for jj = 1:max(size(tau))
    preds = simple_dynamic_link(dp_24, dt, tau(jj));
    [r_temp, p_temp] = corrcoef(preds, do_24); % get the r^2 value
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

% r^2 series
figure
titleout = 'Oith driven by para with month detrend; r^2 values';
plot(tau, r_, '-o', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title(titleout, 'FontSize', 14)
xlabel('time lag', 'FontSize', 12)
ylabel('r^2', 'FontSize', 12)

% predicted and actual
figure
titleout2 = 'Oith driven by para with month detrend; best prediction from AR-1';
plot(xx, do_24, xx, pred_best, 'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', xx)
axis tight
datetick('x', 'mm-dd', 'keeplimits')
grid on
title(titleout2, 'FontSize', 14)
legend('Measured', 'Predicted')
xlabel('Time', 'FontSize', 12)
ylabel('Number of ROIs', 'FontSize', 12)
%%
month = movmean(oith_scaled, 336, 'omitnan');
figure; plot(month)

test = oith_scaled - month;

figure; plot(test)

[ps, freq] = periodogram(test, [1,2]);
figure; loglog(freq, ps)