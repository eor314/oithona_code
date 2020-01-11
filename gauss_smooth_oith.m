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

%%%%% Low pass filter %%%%%
f_cut = 1/(12*24); % cut off frequency as 1/# hours

% FIR filter with hamming window
% d = designfilt('lowpassfir', 'FilterOrder', 10, ...
% 'CutoffFrequency', f_cut, 'DesignMethod', 'window');
% 
% f_delay = (length(d.Coefficients)-1)/2; % compute phase delay

% remove signals by subtracting gaussian filtered data
%gg.w1 = gausswin(24*28); % lunar signal
%gg.w2 = gausswin(24*14);
%gg.w3 = gausswin(24*7);
gg.w4 = gausswin(24);
gg.w4 = gausswin(12);

g_names = fieldnames(gg);
r_names = fieldnames(raw);

detrend.oith = (raw.oith - mean(raw.oith))./std(raw.oith);
detrend.egg = (raw.egg - mean(raw.egg))./std(raw.egg);
detrend.para = (raw.para - mean(raw.para))./std(raw.para);

for ii = 1:length(g_names)
    for jj = 1:length(r_names)
        
        temp = filter(gg.(g_names{ii}), 1, raw.(r_names{jj}));
        temp = (temp - mean(temp))./std(temp);
        detrend.(r_names{jj}) = detrend.(r_names{jj}) - temp;

    end
end

% filter the data
% filt.oith = filter(d, raw.oith);
% filt.para = filter(d, raw.para);
% filt.egg = filter(d, raw.egg);

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
    
    xx = detrend.(names{pairs(ii, 1)});
    yy = detrend.(names{pairs(ii, 2)});
    
    figure; 
    plot(xx, yy, '.', 'MarkerSize', 8)
    grid on
    xlabel(names{pairs(ii,1)}, 'FontSize', 12)
    ylabel(names{pairs(ii,2)}, 'FontSize', 12)
%     title([names{combos(ii,1)}, ' vs. ', names{combos(ii,2)}, ...
%         ' over 2.2 spice'], 'FontSize', 14)
    title(sprintf('%s vs %s, lowpass cutoff %d', ...
        names{pairs(ii, 2)}, names{pairs(ii,1)}, f_cut), 'FontSize', 14)

    %%%%% Cross correlation of eggs and paradinium %%%%%
    % compute the cross correlation keeping xx stationary and lagging
    % yy
    [C, lags] = xcorr(xx, yy, 'coeff');

    figure;
    plot(lags, C, 'LineWidth', 2)
    grid on
    title(sprintf('Cross correlation of %s lagging %s, lowpass cutoff %d', ...
        names{pairs(ii, 2)}, names{pairs(ii,1)}, f_cut), 'FontSize', 14)
%     title(['Cross correlation of ', names{pairs(ii,2)}, ' lagging ', ...
%         names{pairs(ii,1)}, ' lowpass cutoff =', ], 'FontSize', 14)
    xlabel('Lags in hours', 'FontSize', 12)
    ylabel('Correlation coefficent', 'FontSize', 12)
end