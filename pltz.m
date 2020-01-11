clear all 
close all

%%%%% Read data %%%%%
data = dlmread('/Users/Orenstein/Desktop/corrected_counts_091817.txt',',');

data = data(26:end-744, :);

%%%%% Make the corrected hourly counts and bin %%%%%
input1 = {'oith', data(1:end-1, 2) - data(1:end-1, 3) + data(1:end-1, 4)};
input2 = {'para', data(1:end-1, 8) - data(1:end-1, 9) + data(1:end-1, 10)};
input3 = {'egg', data(1:end-1, 5) - data(1:end-1, 6) + data(1:end-1, 7)};
time_frames = [1, 6, 24]; % time frames of interest

mm = 'run_mean'; % define method for averaging
%mm = 'block';

dd = bin_pack(time_frames, mm, input1, input2, input3);

%%%%% Plot the eggs v parasite relative to total oithona %%%%%
oith_tot = dd.hour24.oith + dd.hour24.para + dd.hour24.egg;
eggs = dd.hour24.egg ./ oith_tot;
para = dd.hour24.para ./ oith_tot;

tt = 1:max(size(para));

%tt = datetime(2015, 3, 11) + calhours(1:max(size(para)));

startDate = datenum('03-11-2015 12:00:00 PM');
endDate = datenum('08-01-2015 12:00:00 AM');
%xx = startDate:1/24:endDate;
xx = linspace(startDate, endDate, max(size(para)));
%xx = startDate + tt;
str = datestr(xx, 'mm-dd');

figure;
plot(xx, para, xx, eggs, 'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', xx)
axis tight

%ax.XTick = xx;
datetick('x', 'mm-dd', 'keeplimits')
grid on
xlabel('Date in 2015', 'FontSize', 12)
ylabel('Proportion of total Oithona', 'FontSize', 12)
legend('Paradinium', 'Ovigerous')
title('Relative number of total Oithona', 'FontSize', 14)
%xtickformat('MM-dd')

plot(xx, dd.hour24.oith, xx, dd.hour24.para, xx, dd.hour24.egg, ...
    'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', xx)
axis tight
datetick('x', 'mm-dd', 'keeplimits')
grid on
xlabel('Date in 2015', 'FontSize', 12)
ylabel('Corrected counts', 'FontSize', 12)
legend('Oithona', 'Paradinium', 'Ovigerous')
title('Corrected counts of classes of interest. Summer 2015', 'FontSize', 14)

figure; 
plot(para, eggs, '.', 'MarkerSize', 8)
grid on
xlabel('paradinium', 'FontSize', 12)
ylabel('eggs', 'FontSize', 12)
title('Paradinium vs. Oithona per hour', 'FontSize', 14)

%%%%% Cross correlation of eggs and paradinium %%%%%
% compute the cross correlation keeping eggs stationary and lagging
% paradinium
[C, lags] = xcorr(eggs, para, 'coeff');

figure;
plot(lags, C, 'LineWidth', 2)
grid on
title('Cross correlation of eggs and paradinium', 'FontSize', 14)
xlabel('Lags in hours', 'FontSize', 12)
ylabel('Correlation coefficent', 'FontSize', 12)

%%%%% Spectral plots at 24 hr intervals %%%%%
comp = {'egg', 'para', 'oith'};

for ii = 1:max(size(comp))
    
    in_scaled = (dd.hour24.(comp{ii}) - mean(dd.hour24.(comp{ii}))) / ...
        std(dd.hour24.(comp{ii}));
    
    [ps, freq] = periodogram(in_scaled, [1,2]);
    
    titlestr = ['Spectra ', comp{ii}, ' at 24 hour intervals'];
    
    figure;
    loglog(freq, ps, 'LineWidth', 2)
    grid on
    title(titlestr, 'FontSize', 14)
    xlabel('Log frequency', 'FontSize', 12)
    ylabel('Log power', 'FontSize', 12)
end