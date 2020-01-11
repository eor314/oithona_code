clear all 
close all

% seperate into different bins by spice. Based on looking at plots of spice
% and classes, try cut off at sp = 2.2 [12/17/17]
%
% Filter (low pass and band pass before doing rest of analysis
%% 
% get data
pt_oith = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

[raw, times] = load_oith(pt_oith, '03-11-2015 12:00:00 PM', '08-01-2015 12:00:00 AM');

% names for plotting
names = {'oith', 'para', 'egg'};

%%
% get physical data and compute GSW spice via TEOS-10
pt_phys = '/Users/Orenstein/Documents/documents/SPC/oithona_project/autoss_8172_d7f5_b983.csv';

[phys] = load_sccoos_data(pt_phys);

s_r = gsw_SR_from_SP(phys(:,3));
c_t = gsw_CT_from_t(s_r, phys(:,2), phys(:, 1));

gsw_spice = gsw_spiciness0(s_r, c_t);

sp = gsw_spice(853:end-853); 

sp_dec = decimate(sp, floor(max(size(phys))/max(size(raw.para))));

%%
%%%%% Low pass filter %%%%%
f_cut = 1/(24*7); % cut off frequency as 1/# hours

% FIR filter with hamming window
d = designfilt('lowpassfir', 'FilterOrder', 10, ...
'CutoffFrequency', f_cut, 'DesignMethod', 'window');

% freqz(d, 512, 1)

f_delay = (length(d.Coefficients)-1)/2; % compute phase delay

% filter
for ii = 1:length(names)
    inpt = raw.(names{ii});
    filtdata.(names{ii}) = abs(hilbert(filter(d, detrend(inpt - mean(inpt)))));
end

%%
%%%%% BANDPASS FILTER %%%%%
% this is still wacky not sure what is wrong [12/28/17]

fs = 1/24; % sampling freqency (24 cycles/day for hourly data)
f0 = (1/(24*14)*(1/24)); % center of passband
fn = fs/2; % nyquist frequency assuming sampling of 1 hour
f0n = f0/fn; % normalize the frequency for the notch
n_width = 0.001; % passband width
con = 1.25; %stopband multiplier

bdesign = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    f0n-(con*n_width), f0n-n_width, f0n+n_width, f0n+(con*n_width), 60, 1, 60, fs);

bpass = design(bdesign, 'equiripple');

% filter
for ii = 1:length(names)
    inpt = raw.(names{ii});
    filtdata.(names{ii}) = abs(hilbert(filter(bpass, detrend(inpt - mean(inpt)))));
end

%%

% put it all together
mat = [sp_dec, filtdata.oith, filtdata.para, filtdata.egg];

% find the indicies over and under the limit
over_inds = find(sp_dec > 2.2);
under_inds = find(sp_dec < 2.2);

%over = mat(over_inds, 2:end);
%under = mat(under_inds, 2:end);

under = mat(1:2504, 2:end); % under 2.2 occurs in early part of time series
over = mat(2504:end, 2:end); % over 2.2 occurs in the later part

% indicies of data to consider [1=oith, 2=para, 3=egg]
combos = [1, 2; 1, 3; 3, 2];

for ii = 1:length(combos)
    x_o = over(:, combos(ii,1));
    y_o = over(:, combos(ii,2));

    x_u = under(:, combos(ii,1));
    y_u = under(:, combos(ii,2));

    [C_u, lag_u] = xcorr(x_u, y_u, 1000, 'coeff');
    [C_o, lag_o] = xcorr(x_o, y_o, 1000, 'coeff');

    figure;
    hold on
    set(gca, 'FontSize', 12)

    subplot(2,2,1)
    plot(lag_u, C_u, 'LineWidth', 2)
    grid on
    title(['Cross correlation of ', names{combos(ii,2)}, ' lagging ', ...
        names{combos(ii,1)}, ' under 2.2 spice'], 'FontSize', 14)
    xlabel('Lags in hours', 'FontSize', 12)
    ylabel('Correlation coefficent', 'FontSize', 12)

    subplot(2,2,2)
    plot(lag_o, C_o, 'LineWidth', 2)
    grid on
    title(['Cross correlation of ', names{combos(ii,2)}, ' lagging ', ...
        names{combos(ii,1)}, ' over 2.2 spice'], 'FontSize', 14)
    xlabel('Lags in hours', 'FontSize', 12)
    ylabel('Correlation coefficent', 'FontSize', 12)

    subplot(2,2,3)
    plot(x_u, y_u, '.', 'MarkerSize', 8)
    grid on
    xlabel(names{combos(ii,1)}, 'FontSize', 12)
    ylabel(names{combos(ii,2)}, 'FontSize', 12)
    title([names{combos(ii,1)}, ' vs. ', names{combos(ii,2)}, ...
        ' under 2.2 spice'], 'FontSize', 14)

    subplot(2,2,4)
    plot(x_o, y_o, '.', 'MarkerSize', 8)
    grid on
    xlabel(names{combos(ii,1)}, 'FontSize', 12)
    ylabel(names{combos(ii,2)}, 'FontSize', 12)
    title([names{combos(ii,1)}, ' vs. ', names{combos(ii,2)}, ...
        ' over 2.2 spice'], 'FontSize', 14)
end