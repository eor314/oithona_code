clear all 
close all

% Play with MATLABs AR toolbox

%% Load the data

pt_oith = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

[raw, tt] = load_oith(pt_oith, '03-11-2015 12:00:00 PM', '08-01-2015 12:00:00 AM');

%% Generate a bunch of different plots at different lags

smooth = movmean(raw.para, 24);

%%

for ii = linspace(24, 7*24, 7)
    
    x1 = smooth(1:end-(ii-1));
    x2 = smooth(ii:end);
    
    figure; plot(x1, x2, '.', 'MarkerSize', 6)
    title(sprintf('Lag = %d', ii), 'FontSize', 14)
end

%% Sample autocorrelation sequence out to 50
[xc, lags] = xcorr(smooth, 100, 'coeff');

figure; 
stem(lags(101:end), xc(101:end), 'filled')
xlabel('Lag')
ylabel('ACF')

%%
[arcoef, E, K] = aryule(raw.egg, 50);
pacf = -K;

figure;
stem(pacf, 'filled')
xlim([1 50])
uconf = 1.96/sqrt(length(raw.para));
lconf = -uconf;
hold on
plot([1 50], [1 1]'*[lconf uconf], 'r')
grid on

%% 
fil = filter(1, arcoef, raw.egg);

figure;
plot(fil)
title('filtered data')

[C, lags] = xcorr(raw.eggs, fil, 1000, 'coeff');

%% 
fil = filter(1, arcoef(1:2), raw.egg);

[C, lags] = xcorr(raw.eggs, fil, 1000, 'coeff');
