clear all 
close all

%% 
% get data
pt_oith = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

[raw, tt] = load_oith(pt_oith, '03-11-2015 12:00:00 PM', '08-01-2015 12:00:00 AM');

smooth_oith = movmean(raw.oith, 24);
smooth_egg = movmean(raw.egg, 24);
smooth_para = movmean(raw.para, 24);

% chunk into lunar cycles (hard coded)
cycles = struct();
cycles(1).oith = smooth_oith(191:886);
cycles(1).egg = smooth_egg(191:886);
cycles(1).para = smooth_para(191:886);

cycles(2).oith = smooth_oith(888:1607);
cycles(2).egg = smooth_egg(888:1607);
cycles(2).para = smooth_para(888:1607);

cycles(3).oith = smooth_oith(1608:2303);
cycles(3).egg = smooth_egg(1608:2303);
cycles(3).para = smooth_para(1608:2303);

cycles(4).oith = smooth_oith(2304:3023);
cycles(4).egg = smooth_egg(2304:3023);
cycles(4).para = smooth_para(2304:3023);

mths = {'mar', 'apr', 'jun', 'jul'};
names = {'oith', 'para', 'egg'};

%% plot xcorrs
combos = [1, 2; 3, 2; 1, 2];

for ii = 1:size(combos,1)
    
    figure;
    for jj = 1:length(mths)

        xx = cycles(jj).(names{combos(ii, 1)});
        yy = cycles(jj).(names{combos(ii, 2)});
        
        [C, lags] = xcorr(xx, yy, 'coeff');
        
        subplot(2,2,jj)
        plot(lags, C, 'LineWidth', 2)
        grid on
        title(['Cross correlation of ', names{combos(ii,2)}, ' lagging ', ...
            names{combos(ii,1)}, ' in ', mths{jj}], 'FontSize', 14)
        xlabel('Lags in hours', 'FontSize', 12)
        ylabel('Correlation coefficent', 'FontSize', 12)
    end
    

end
        