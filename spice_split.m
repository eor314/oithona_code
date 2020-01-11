 clear all 
close all

% seperate into different bins by spice. Based on looking at plots of spice
% and classes, try cut off at sp = 2.2 [12/17/17]
%% 
% get data
pt_oith = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

[raw, times] = load_oith(pt_oith, '03-11-2015 12:00:00 PM', '08-01-2015 12:00:00 AM');

%%
% get physical data and compute GSW spice via TEOS-10
pt_phys = '/Users/Orenstein/Documents/documents/SPC/oithona_project/autoss_8172_d7f5_b983.csv';

[phys] = load_sccoos_data(pt_phys);

s_r = gsw_SR_from_SP(phys(:,3));
c_t = gsw_CT_from_t(s_r, phys(:,2), phys(:, 1));

gsw_spice = gsw_spiciness0(s_r, c_t);

sp = gsw_spice(853:end-853); 

sp_dec = decimate(sp, floor(max(size(phys))/max(size(raw.para))));

% white.oith = (raw.oith - mean(raw.oith));
% white.para = (raw.para - mean(raw.para));
% white.egg = (raw.egg - mean(raw.egg));

% smoothed stuff to 24 hrs
% bin_para = movmean(raw.para, 24, 'omitnan');
% bin_egg = movmean(raw.egg, 24, 'omitnan');
% bin_oith = movmean(raw.oith, 24, 'omitnan');

% remove mean and hilbert transform
bin_para = abs(hilbert(detrend(raw.para - mean(raw.para))));
bin_egg = abs(hilbert(detrend(raw.egg - mean(raw.egg))));
bin_oith = abs(hilbert(detrend(raw.oith - mean(raw.oith))));

% whiten
% bin_para = (bin_para - mean(bin_para))./std(bin_para);
% bin_egg = (bin_egg - mean(bin_egg))./std(bin_egg);
% bin_oith = (bin_oith - mean(bin_oith))./std(bin_oith);

% put it all together
mat = [sp_dec, bin_oith, bin_para, bin_egg];

% find the indicies over and under the limit
over_inds = find(sp_dec > 2.2);
under_inds = find(sp_dec < 2.2);

%over = mat(over_inds, 2:end);
%under = mat(under_inds, 2:end);

under = mat(1:2504, 2:end); % under 2.2 occurs in early part of time series
over = mat(2504:end, 2:end); % over 2.2 occurs in the later part

% names for plotting
names = {'oith', 'para', 'egg'};

% indicies of data to consider [1=oith, 2=para, 3=egg]
%combos = [1, 2; 2, 3; 1, 3; 3, 2];

%combos = [2, 3; 3, 2]; % para and eggs
%combos = [1,2; 2,1]; % oith and para
combos = [1, 3; 3, 1]; % oith and eggs
%% 

for ind = 1:length(combos)

    xx = over(:, combos(ind, 1));
    yy = over(:, combos(ind, 2));
    
    figure; 
    plot(xx, yy, '.', 'MarkerSize', 8)
    grid on
    xlabel(names{combos(ind,1)}, 'FontSize', 12)
    ylabel(names{combos(ind,2)}, 'FontSize', 12)
    title([names{combos(ind,1)}, ' vs. ', names{combos(ind,2)}, ...
        ' over 2.2 spice'], 'FontSize', 14)

    %%%%% Cross correlation of eggs and paradinium %%%%%
    % compute the cross correlation keeping xx stationary and lagging
    % yy
    [C, lags] = xcorr(xx, yy, 1000, 'coeff');

    figure;
    plot(lags, C, 'LineWidth', 2)
    grid on
    title(['Cross correlation of ', names{combos(ind,2)}, ' lagging ', ...
        names{combos(ind,1)}, ' over 2.2 spice'], 'FontSize', 14)
    xlabel('Lags in hours', 'FontSize', 12)
    ylabel('Correlation coefficent', 'FontSize', 12)
end
%% 

for ind = 1:length(combos)

    xx = under(:, combos(ind, 1));
    yy = under(:, combos(ind, 2));
    
    figure; 
    plot(xx, yy, '.', 'MarkerSize', 8)
    grid on
    xlabel(names{combos(ind,1)}, 'FontSize', 12)
    ylabel(names{combos(ind,2)}, 'FontSize', 12)
    title([names{combos(ind,1)}, ' vs. ', names{combos(ind,2)}, ...
        ' under 2.2 spice'], 'FontSize', 14)

    %%%%% Cross correlation of eggs and paradinium %%%%%
    % compute the cross correlation keeping xx stationary and lagging
    % yy
    [C, lags] = xcorr(xx, yy, 1000, 'coeff');

    figure;
    plot(lags, C, 'LineWidth', 2)
    grid on
    title(['Cross correlation of ', names{combos(ind,2)}, ' lagging ', ...
        names{combos(ind,1)}, ' under 2.2 spice'], 'FontSize', 14)
    xlabel('Lags in hours', 'FontSize', 12)
    ylabel('Correlation coefficent', 'FontSize', 12)
end

%% put stuff on same plot
combo = [2, 3];
   
x_o = over(:, combo(1));
y_o = over(:, combo(2));

x_u = under(:, combo(1));
y_u = under(:, combo(2));

[C_u, lag_u] = xcorr(x_u, y_u, 1000, 'coeff');
[C_o, lag_o] = xcorr(x_o, y_o, 1000, 'coeff');

figure;
hold on
set(gca, 'FontSize', 12)

subplot(2,2,1)
plot(lag_u, C_u, 'LineWidth', 2)
grid on
title(['Cross correlation of ', names{combo(2)}, ' lagging ', ...
    names{combo(1)}, ' under 2.2 spice'], 'FontSize', 14)
xlabel('Lags in hours', 'FontSize', 12)
ylabel('Correlation coefficent', 'FontSize', 12)

subplot(2,2,2)
plot(lag_o, C_o, 'LineWidth', 2)
grid on
title(['Cross correlation of ', names{combo(2)}, ' lagging ', ...
    names{combo(1)}, ' over 2.2 spice'], 'FontSize', 14)
xlabel('Lags in hours', 'FontSize', 12)
ylabel('Correlation coefficent', 'FontSize', 12)

subplot(2,2,3)
plot(x_u, y_u, '.', 'MarkerSize', 8)
grid on
xlabel(names{combo(1)}, 'FontSize', 12)
ylabel(names{combo(2)}, 'FontSize', 12)
title([names{combo(1)}, ' vs. ', names{combo(2)}, ...
    ' under 2.2 spice'], 'FontSize', 14)

subplot(2,2,4)
plot(x_o, y_o, '.', 'MarkerSize', 8)
grid on
xlabel(names{combo(1)}, 'FontSize', 12)
ylabel(names{combo(2)}, 'FontSize', 12)
title([names{combo(1)}, ' vs. ', names{combo(2)}, ...
    ' over 2.2 spice'], 'FontSize', 14)

%%
%%%%% Run AR-1 model %%%%%

combos = [3, 1; 2, 3; 2, 1];
% define lags and pre-allocate space for r2 and significance metrics  
dt = 1;
tau = 1:(dt/4):(15/dt);

for ii = 1:length(combos)
    
    x_o = over(:, combos(ii, 1));
    y_o = over(:, combos(ii, 2));
    
    x_u = under(:, combos(ii, 1));
    y_u = under(:, combos(ii, 2));
    
    r_o = zeros(size(tau));
    p_u = zeros(size(tau));
    
    r_u = zeros(size(tau));
    p_o = zeros(size(tau));

    %%%%% do the prediction for over %%%%%
    for jj = 1:max(size(tau))
        preds = simple_dynamic_link(x_o, dt, tau(jj));
        [r_temp, p_temp] = corrcoef(preds, y_o); % get the r^2 value
        r_o(jj) = r_temp(1,2);
        p_o(jj) = p_temp(1,2);
        if jj > 1
            if r_o(jj) > r_o(jj-1)
                pred_best = preds;
            end
        else
            pred_best = preds;
        end
        clear temp
    end
    
    % predicted and actual over 2.2
    figure
    titleout2 = [names{combos(ii,2)}, ' driven by ', names{combos(ii,1)}, ...
        ' over 2.2: actual data and best AR-1 prediction'];
    plot(times.xx(2504:end), y_o, times.xx(2504:end), pred_best, 'LineWidth', 2)
    labs = get(gca,'XLim');
    set(gca, 'FontSize', 12)
    set(gca, 'XTick', times.xx(2504:end))
    axis tight
    datetick('x', 'mm-dd', 'keeplimits')
    grid on
    title(titleout2, 'FontSize', 14)
    legend('Measured', 'Predicted')
    xlabel('Time', 'FontSize', 12)
    ylabel('Scaled Number of ROIs', 'FontSize', 12)

    % r^2 series over 2.2
    figure
    titleout = [names{combos(ii,2)}, ' driven by ', names{combos(ii,1)}, ...
        ' over 2.2: r^2 values'];
    plot(tau, r_o, '-o', 'LineWidth', 2, 'MarkerSize', 8)
    set(gca, 'FontSize', 12)
    grid on
    title(titleout, 'FontSize', 14)
    xlabel('time lag', 'FontSize', 12)
    ylabel('r^2', 'FontSize', 12)

    % Print the maximum r and associated tau over 2.2
    [xx, ind] = max(r_o);

    sprintf([names{combos(ii,2)}, ' driven by ', names{combos(ii,1)}, '\n',...
        'r^2 max over 2.2 = %0.5g \n tau = %0.5g \n'], xx, tau(ind))
    
    %%%%% do the prediction for under %%%%%
    for jj = 1:max(size(tau))
        preds = simple_dynamic_link(x_u, dt, tau(jj));
        [r_temp, p_temp] = corrcoef(preds, y_u); % get the r^2 value
        r_u(jj) = r_temp(1,2);
        p_u(jj) = p_temp(1,2);
        if jj > 1
            if r_u(jj) > r_u(jj-1)
                pred_best = preds;
            end
        else
            pred_best = preds;
        end
        clear temp
    end

    % predicted and actual under 2.2
    figure
    titleout2 = [names{combos(ii,2)}, ' driven by ', names{combos(ii,1)}, ...
        ' under 2.2: actual data and best AR-1 prediction'];
    plot(times.xx(1:2504), y_u, times.xx(1:2504), pred_best, 'LineWidth', 2)
    labs = get(gca,'XLim');
    set(gca, 'FontSize', 12)
    set(gca, 'XTick', times.xx(1:2504))
    axis tight
    datetick('x', 'mm-dd', 'keeplimits')
    grid on
    title(titleout2, 'FontSize', 14)
    legend('Measured', 'Predicted')
    xlabel('Time', 'FontSize', 12)
    ylabel('Scaled Number of ROIs', 'FontSize', 12)

    % r^2 series over 2.2
    figure
    titleout = [names{combos(ii,2)}, ' driven by ', names{combos(ii,1)}, ...
        ' under 2.2: r^2 values'];
    plot(tau, r_u, '-o', 'LineWidth', 2, 'MarkerSize', 8)
    set(gca, 'FontSize', 12)
    grid on
    title(titleout, 'FontSize', 14)
    xlabel('time lag', 'FontSize', 12)
    ylabel('r^2', 'FontSize', 12)

    % Print the maximum r and associated tau under 2.2 
    [xx, ind] = max(r_u);

    sprintf([names{combos(ii,2)}, ' driven by ', names{combos(ii,1)}, '\n',...
        'r^2 max under 2.2 = %0.5g \n tau = %0.5g \n'], xx, tau(ind))
end