clear all 
close all

% lowpass filter the data at a cutoff frequency and look at the
% xcorrelations and data plotted against one another

%% 
% get data
pt_oith = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

[raw, tt] = load_oith(pt_oith, '03-11-2015 12:00:00 PM', '08-01-2015 12:00:00 AM');

smooth.oith = movmean(raw.oith, 24);
smooth.egg = movmean(raw.egg, 24);
smooth.para = movmean(raw.para, 24);

%% 
%%%%% Notch parameters %%%%%
f0 = 1/(24*14); % notch frequency 1
f1 = 1/12; % notch frequency 2
fn = 1/2; % nyquist frequency assuming sampling of 1 hour
f0n = f0/fn; % normalize the frequency for the notch
f1n = f1/fn;
n_width = 0.001; % notch width

%%
%%%%% Run notchj filter %%%%%
dd = fdesign.notch(8, f0n, 10);
Hd = design(dd, 'Systemobject', true);

% filter the data
filt.oith = Hd(smooth.oith);
filt.para = Hd(smooth.para);
filt.egg = Hd(smooth.egg);

% filter smoothed data
% filt.oith = filter(d, smooth_oith);
% filt.para = filter(d, smooth_para);
% filt.egg = filter(d, smooth_egg);

% plot it
%{
figure; plot(tt.xx - f_delay/24, filt.oith, tt.xx, raw.oith)
title('oith')
figure; plot(tt.xx - f_delay/24, filt.para, tt.xx, raw.para)
title('para')
figure; plot(tt.xx - f_delay/24, filt.egg, tt.xx, raw.egg)
title('egg')
%}
names = {'oith', 'para', 'egg'};
%%
pairs = [1, 2; 3, 2; 1, 3];

out = struct();

for ii = 1:3
    
    xx = filt.(names{pairs(ii, 1)});
    yy = filt.(names{pairs(ii, 2)});
    
    figure; 
    plot(xx, yy, '.', 'MarkerSize', 8)
    grid on
    xlabel(names{pairs(ii,1)}, 'FontSize', 12)
    ylabel(names{pairs(ii,2)}, 'FontSize', 12)
%     title([names{combos(ii,1)}, ' vs. ', names{combos(ii,2)}, ...
%         ' over 2.2 spice'], 'FontSize', 14)
    title(sprintf('%s vs %s, notch cutoff %d', ...
        names{pairs(ii, 2)}, names{pairs(ii,1)}, f0n), 'FontSize', 14)

    %%%%% Cross correlation of eggs and paradinium %%%%%
    % compute the cross correlation keeping xx stationary and lagging
    % yy
    [C, lags] = xcorr(xx, yy, 'coeff');

    figure;
    plot(lags, C, 'LineWidth', 2)
    grid on
    title(sprintf('Cross correlation of %s lagging %s, notch cutoff %d', ...
        names{pairs(ii, 2)}, names{pairs(ii,1)}, f0n), 'FontSize', 14)
%     title(['Cross correlation of ', names{pairs(ii,2)}, ' lagging ', ...
%         names{pairs(ii,1)}, ' lowpass cutoff =', ], 'FontSize', 14)
    xlabel('Lags in hours', 'FontSize', 12)
    ylabel('Correlation coefficent', 'FontSize', 12)
end