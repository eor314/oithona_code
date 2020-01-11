clear all 
close all

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

% what to plot
str = 'para';

% smoothed paradinium
bin = movmean(raw.(str), 24, 'omitnan');

%% plot spice v. para
% remove some data from either side to make series resample nicely to same
% length as oithona data [hardcoded]
sp = gsw_spice(853:end-853); 

sp_dec = decimate(sp, floor(max(size(phys))/max(size(raw.(str)))));

figure; 
plot(sp_dec, bin, '.', 'MarkerSize', 12)
set(gca, 'FontSize', 12)
title(['spice v. 24 hr smoothed ', str], 'FontSize', 14)
xlabel('Spice', 'FontSize', 12)
ylabel(['Number ', str], 'FontSize', 12)
grid on

%% plot temp v. para
temp = phys(853:end-853, 2); 

temp_dec = decimate(temp, floor(max(size(phys))/max(size(raw.(str)))));

figure; 
plot(temp_dec, bin, '.', 'MarkerSize', 12)
set(gca, 'FontSize', 12)
title(['temp v. 24 hr smoothed ', str], 'FontSize', 14)
xlabel('temp', 'FontSize', 12)
ylabel(['Number ', str], 'FontSize', 12)
grid on

%% plot sal v. para
sal = phys(853:end-853, 3); 

sal_dec = decimate(sal, floor(max(size(phys))/max(size(raw.(str)))));
%sal_dec = downsample(sal, floor(max(size(phys))/max(size(raw.(str)))));

figure; 
plot(sal_dec, bin, '.', 'MarkerSize', 12)
set(gca, 'FontSize', 12)
title(['sal v. 24 hr smoothed ', str], 'FontSize', 14)
xlabel('sal', 'FontSize', 12)
ylabel(['Number ', str], 'FontSize', 12)
grid on

%% plot pressure v. para
press = phys(853:end-853, 1); 

press_dec = decimate(press, floor(max(size(phys))/max(size(raw.(str)))));

figure; 
plot(press_dec, bin, '.', 'MarkerSize', 8)
set(gca, 'FontSize', 12)
title(['press v. 24 hr smoothed ', str], 'FontSize', 14)
xlabel('press', 'FontSize', 12)
ylabel(['Number ', str], 'FontSize', 12)
grid on
