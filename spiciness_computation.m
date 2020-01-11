clear all
%% read data and smooth
ff = readtable('/Users/Orenstein/Documents/documents/SPC/oithona_project/autoss_8172_d7f5_b983.csv');

ff = ff(1:34345, :); % just select till august

% smooth over 24 hour periods. Data collected at 6 minute intervals
ff.salinity_1e_3_ = movmean(ff.salinity_1e_3_, 240, 'omitnan');
ff.temperature_celsius_ = movmean(ff.temperature_celsius_, 240, 'omitnan');
ff.pressure_dbar_ = movmean(ff.pressure_dbar_, 240, 'omitnan');

%% get rid of any bad data points according to data QC flag
p_flag = ff.pressure_flagPrimary==1;
s_flag = ff.salinity_flagPrimary==1;
t_flag = ff.temperature_flagPrimary==1;

flags = p_flag + s_flag + t_flag;
flags = (flags == 3);

xx = ff(flags, :);

%% get rid of any questionable salinity points outside of 2 STDs
sd_s = std(xx.salinity_1e_3_);
mu_s = mean(xx.salinity_1e_3_);

under = xx.salinity_1e_3_ > (mu_s - 2*sd_s);
over = xx.salinity_1e_3_ < (mu_s + 2*sd_s);

ou = under + over;
ou = (ou == 2);

xx = xx(ou, :);

clear -VARS under over ou

%% get rid of any questionable temperature points outside of 2 STDs
sd_t = std(xx.temperature_celsius_);
mu_t = mean(xx.temperature_celsius_);

under = xx.temperature_celsius_ > (mu_t - 2*sd_t);
over = xx.temperature_celsius_ < (mu_t + 2*sd_t);

ou = under + over;
ou = (ou == 2);

xx = xx(ou, :);

clear -VARS under over ou

%% get rid of any questionable pressure points outside of 2 STDs
sd_p = std(xx.pressure_dbar_);
mu_p = mean(xx.pressure_dbar_);

under = xx.pressure_dbar_ > (mu_p - 2*sd_p);
over = xx.pressure_dbar_ < (mu_p + 2*sd_p);

ou = under + over;
ou = (ou == 2);

xx = xx(ou, :);

clear -VARS under over ou

times = xx.time_UTC_;
out = xx{:, {'pressure_dbar_', 'temperature_celsius_', 'salinity_1e_3_'}};

%% get datetime for plotting
startDate = datenum('03-11-2015 12:00:00 PM');
endDate = datenum('08-01-2015 12:00:00 AM');
dd = linspace(startDate, endDate, max(size(xx)));
str = datestr(dd, 'mm-dd');

%% try flament spice
flam_sp = flament_spice(out(:, 2), out(:, 3));

figure;
plot(dd, flam_sp)
labs = get(gca,'XLim');
set(gca, 'XTick', dd)
axis tight
datetick('x', 'mm-dd', 'keeplimits')
grid on
title('Flament Spice from PSU and temp, 24 hr smooth')
xlabel('Time')
ylabel('Spiciness')

%% GSW spice via TEOS-10
s_r = gsw_SR_from_SP(out(:,3));
c_t = gsw_CT_from_t(s_r, out(:,2), out(:, 1));

gsw_spice = gsw_spiciness0(s_r, c_t);

figure;
plot(dd, gsw_spice)
labs = get(gca,'XLim');
set(gca, 'XTick', dd)
axis tight
datetick('x', 'mm-dd', 'keeplimits')
grid on
title('GSW Spice from ref. sal. and consv. temp., 24 hr smooth')
xlabel('Time')
ylabel('Spiciness')

%% Get oithona data and process
ptf = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

data = dlmread(ptf, ',');
nz_front = find(data(:, 2), 1, 'first');
nz_end = find(data(:, 2), 1, 'last');

% pull out the relevent data
data = data(nz_front:nz_end, :); 
oith = data(1:end, 2) - data(1:end, 3) + data(1:end, 4);
para = data(1:end, 8) - data(1:end, 9) + data(1:end, 10);
egg = data(1:end, 5) - data(1:end, 6) + data(1:end, 7);

oith_tot = oith + para + egg;
para_rel = para ./ oith_tot;
egg_rel = egg ./ oith_tot;

% smooth and scale each time series
oith_tot = movmean(oith_tot, 24, 'omitnan');
%oith_tot = (oith_tot - mean(oith_tot)) / std(oith_tot);
para_rel = movmean(para_rel, 24, 'omitnan');
egg_rel = movmean(egg_rel, 24, 'omitnan');

oith = movmean(oith, 24, 'omitnan');
oith = (oith - mean(oith)) / std(oith);

para = movmean(para, 24, 'omitnan');
para = (para - mean(para)) / std(para);

egg = movmean(egg, 24, 'omitnan');
egg = (egg - mean(egg)) / std(egg);

%% plot everything with spiciness
dd2 = linspace(startDate, endDate, max(size(oith)));

figure;
plot(dd, gsw_spice, 'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', dd)
datetick('x', 'mm-dd', 'keeplimits')
grid on
xlabel('Date in 2015', 'FontSize', 12)
ylabel('Spiciness', 'FontSize', 12)
hold on

yyaxis right
plot(dd2, oith, '-.r', 'LineWidth', 2)
plot(dd2, para, '--g', 'LineWidth', 2)
plot(dd2, egg, ':k', 'LineWidth', 2)
ylabel('# ROI', 'FontSize', 12)

ticklabs = get(gca, 'YTickLabels');
ticklabs_new = cell(size(ticklabs));
for ii = 1:length(ticklabs)
    ticklabs_new{ii} = ['\color{black} ' ticklabs{ii}];
end

set(gca, 'YTickLabel', ticklabs_new)
set(gca, 'YColor', 'black')

title('GSW spiciness, Oithona from SPC, 1-day bin', 'FontSize', 14)
legend('Spiciness', 'Oithona', 'Parasitized', 'Ovigerious')

%% plot total oithona with spiciness
figure;
plot(dd, gsw_spice, 'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', dd)
datetick('x', 'mm-dd', 'keeplimits')
grid on
xlabel('Date in 2015', 'FontSize', 12)
ylabel('Spiciness', 'FontSize', 12)
set(gca, 'FontSize', 12)
hold on

yyaxis right
plot(dd2, oith_tot, '-.r', 'LineWidth', 2)
ylabel('# ROI', 'FontSize', 12)

ticklabs = get(gca, 'YTickLabels');
ticklabs_new = cell(size(ticklabs));
for ii = 1:length(ticklabs)
    ticklabs_new{ii} = ['\color{black} ' ticklabs{ii}];
end

set(gca, 'YTickLabel', ticklabs_new)
set(gca, 'YColor', 'black')
set(gca, 'FontSize', 12)

title('GSW spiciness, Total Oithona, 1-day bin', 'FontSize', 14)
legend({'Spiciness', 'Oithona'}, 'Location', 'northwest', 'FontSize', 12)
axis tight

%% plot relative eggs and paradinium with spiciness
figure;
plot(dd, gsw_spice, 'LineWidth', 2)
labs = get(gca,'XLim');
set(gca, 'XTick', dd)
datetick('x', 'mm-dd', 'keeplimits')
grid on
xlabel('Date in 2015', 'FontSize', 12)
ylabel('Spiciness', 'FontSize', 12)
hold on

yyaxis right
plot(dd2, para_rel, '--g', 'LineWidth', 2)
plot(dd2, egg_rel, ':r', 'LineWidth', 2)
ylabel('Percent', 'FontSize', 12, 'Color', 'k')

ticklabs = get(gca, 'YTickLabels');
ticklabs_new = cell(size(ticklabs));
for ii = 1:length(ticklabs)
    ticklabs_new{ii} = ['\color{black} ' ticklabs{ii}];
end

set(gca, 'YTickLabel', ticklabs_new)
set(gca, 'YColor', 'black')

title('GSW spiciness, % parasitize and ovigerous, 1-day bin', 'FontSize', 14)
legend('Spiciness', '% Parasitized', '% Ovigerious')
axis tight

%% interpolate oithona data and do xcorr and regression
o_interp = interp1(dd2, oith, dd); % linear interpolation

[C, lags] = xcorr(o_interp', gsw_spice, 'coeff');

figure;
plot(lags, C, 'LineWidth', 2)
grid on
title('Cross correlation of total oithona and spiciness', 'FontSize', 14)
xlabel('Lags', 'FontSize', 12)
ylabel('Correlation coefficent', 'FontSize', 12)
set(gca, 'FontSize', 12)
