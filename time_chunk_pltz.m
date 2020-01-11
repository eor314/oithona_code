clear all 
close all

%% 
% get data
pt_oith = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

[raw, tt] = load_oith(pt_oith, '03-11-2015 12:00:00 PM', '08-01-2015 12:00:00 AM');

smooth_oith = movmean(raw.oith, 24);
smooth_egg = movmean(raw.egg, 24);
smooth_para = movmean(raw.para, 24);

% trim a day off either end to get exact weeks
smooth_oith = smooth_oith(25:end-24);
smooth_egg = smooth_egg(25:end-24);
smooth_para = smooth_para(25:end-24);

% time period (hours x days)
amt = 24*14;

% chunk the data
bin.oith = reshape(smooth_oith, amt, max(size(smooth_oith))/(amt));
bin.egg = reshape(smooth_egg, amt, max(size(smooth_egg))/(amt));
bin.para = reshape(smooth_para, amt, max(size(smooth_para))/(amt));

names = {'oith', 'para', 'egg'};
figure; plot(bin.('oith'), bin.('para'), '.', 'MarkerSize', 12)

%% compute all xcorr and regressions
pairs = [2, 1; 2, 3; 3, 1];

out = struct();

for jj = 1:3
    
    xx = bin.(names{pairs(jj, 1)});
    yy =bin.(names{pairs(jj, 2)});
    
    r_ = zeros(max(size(smooth_para))/amt, 2);
    slope = zeros(max(size(smooth_para))/amt, 1);
    
    for ii = 1:max(size(smooth_para))/amt
        xx_seg = xx(:, ii);
        yy_seg = yy(:, ii);
        
        % compute xcorr
        [C, lags] = xcorr(yy_seg, xx_seg, 'coeff');
        [mx, ind] = max(C);
        
        r_(ii,1) = mx;
        r_(ii,2) = lags(ind);
        
        % linear regresssion
        pp = polyfit(xx_seg, yy_seg, 1);
        slope(ii) = pp(1);
    end
    
    str = [names{pairs(jj, 1)},' lagging ', names{pairs(jj, 2)}];
    out(jj).name = str;
    out(jj).xcorr = r_;
    out(jj).slope = slope;
end

%% plot all for one pair
pair = [2, 3];

xx = bin.(names{pair(1)});
yy = bin.(names{pair(2)});

for ii = 1:max(size(smooth_para))/amt
    xx_seg = xx(:, ii);
    yy_seg = yy(:, ii);

    % compute xcorr
    [C, lags] = xcorr(yy_seg, xx_seg, 'coeff');
    
    figure; plot(lags, C)
    title(sprintf('time chunk %d', ii))
end

