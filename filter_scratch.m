
clear all 
close all

%%%%% Input parameters %%%%%

% where to find comma seperated data
ptf = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

data = dlmread(ptf, ',');

% get first and last nonzero indicies
nz_front = find(data(:, 2), 1, 'first');
nz_end = find(data(:, 2), 1, 'last');

% pull out the relevent data
data = data(nz_front:nz_end, :); % hard coded to truncate zeros
raw.oith = data(1:end, 2) - data(1:end, 3) + data(1:end, 4);
raw.para = data(1:end, 8) - data(1:end, 9) + data(1:end, 10);
raw.egg = data(1:end, 5) - data(1:end, 6) + data(1:end, 7);

% names
names = {'oith', 'para', 'egg'};

%%%%% Date and time labels for plotting %%%%%
startDate = datenum('03-11-2015 12:00:00 PM');
endDate = datenum('08-01-2015 12:00:00 AM');
xx = linspace(startDate, endDate, max(size(data)));
str = datestr(xx, 'mm-dd');

%%
%%%%% Power Spectra %%%%%
% sampling frequency of 1 hour

for ii = 1:3
    in_ser = raw.(names{ii});
    [pxx, freq] = pwelch(in_ser-mean(in_ser), 1024, [], [], 1);
    figure;
    semilogx(freq, 10*log10(pxx))
    title(['Welch PSD ', names{ii}, ' win = 128, 50% overlap'], 'FontSize', 14)
    xlabel('Frequency (1/hours)', 'FontSize', 12)
    ylabel('Power', 'FontSize', 12)
    set(gca, 'FontSize', 12)
    grid on
end

%%
%%%%% Low pass filter %%%%%
f_cut = 1/(12*24); % cut off frequency as 1/# hours

% FIR filter with hamming window
d = designfilt('lowpassfir', 'FilterOrder', 10, ...
'CutoffFrequency', f_cut, 'DesignMethod', 'window');

% FIR filter with equiripple
% d = designfilt('lowpassfir',  ...
% 'CutoffFrequency', f_cut, 'DesignMethod', 'equiripple');

f_delay = (length(d.Coefficients)-1)/2; % compute phase delay

% filter the data
filt.oith = filter(d, raw.oith);
filt.para = filter(d, raw.para);
filt.egg = filter(d, raw.egg);

% plot it
figure; plot(xx - f_delay/24, filt.oith, xx, raw.oith)
title('oith')
figure; plot(xx - f_delay/24, filt.para, xx, raw.para)
title('para')
figure; plot(xx - f_delay/24, filt.egg, xx, raw.egg)
title('egg')

%% 
%%%%% Notch parameters %%%%%
f0 = 1/24; % notch frequency 1
f1 = 1/12; % notch frequency 2
fn = 1/2; % nyquist frequency assuming sampling of 1 hour
f0n = f0/fn; % normalize the frequency for the notch
f1n = f1/fn;
n_width = 0.001; % notch width

%%
%%%%% Butterworth band stop %%%%%
wg1 = [f0n - n_width, f0n + n_width];
[b1, a1] = butter(2, wg1, 'stop');

wg2 = [f1n - n_width, f1n + n_width];
[b2, a2] = butter(2, wg2, 'stop');

a0 = conv(a1, a2);
b0 = conv(b1, b2);

test = filter(b0, a0, raw.oith);

[pxx, freq] = pwelch(test - mean(test), 128, [], [], 1);
figure; semilogx(freq, 10*log10(pxx))
grid on 
title('Butterworth band stop')

%%
%%%%% Notch using fdesign %%%%%
dd = fdesign.notch(8, f0n, 10);
Hd = design(dd, 'Systemobject', true);
test = Hd(raw.oith);

[pxx, freq] = pwelch(test - mean(test), 128, [], [], 1);
figure; semilogx(freq, 10*log10(pxx))
grid on
title('fdesign notch')