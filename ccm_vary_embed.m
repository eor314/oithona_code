clear all
close all

%% Load data and get parameters
pt_oith = '/media/storage/learning_files/oithona_project/corrected_counts_091817.txt';

[raw, times] = load_oith(pt_oith, '03-11-2015 12:00:00 PM', '08-01-2015 12:00:00 AM');

%% Low pass filter the time series
%%%%% Low pass filter %%%%%
f_cut = 1/(24*7); % cut off frequency as 1/# hours

% FIR filter with hamming window
d = designfilt('lowpassfir', 'FilterOrder', 10, ...
'CutoffFrequency', f_cut, 'DesignMethod', 'window');

f_delay = (length(d.Coefficients)-1)/2; % compute phase delay

% hilbert transform
filt.oith = abs(hilbert(filter(d, detrend((raw.oith - mean(raw.oith))./std(raw.oith)))));
filt.para = abs(hilbert(filter(d, detrend((raw.para - mean(raw.para))./std(raw.para)))));
filt.egg = abs(hilbert(filter(d, detrend((raw.egg - mean(raw.egg))./std(raw.egg)))));


test.oith = abs(hilbert(detrend((raw.oith - mean(raw.oith))./std(raw.oith))));
test.para = abs(hilbert(detrend((raw.para - mean(raw.para))./std(raw.para))));
test.egg = abs(hilbert(detrend((raw.egg - mean(raw.egg))./std(raw.egg))));

%% try the CCM varying just embedding dimension
%[corrs , predY, predX, origY, origX] = sugi_CCM(filt.para, filt.oith, tau , E);
names = {'oith', 'egg', 'para'};
pairs = [1, 2; 1, 3; 2, 3];

all_out = struct();
tau = 5;

for ii = 1:length(pairs)
    out = struct();
    temp_corr = zeros(2, 20);
    for ee = 1:20
        [corrs, ~] = sugi_CCM(test.(names{pairs(ii,1)}), test.(names{pairs(ii,2)}), ...
            tau, ee);
        temp_corr(:, ee) = corrs;
    end

    titlestr1 = sprintf([names{pairs(ii,1)}, ' and ', ...
        names{pairs(ii,2)}, ' E = 1:%d, %s = %d'], ee, '$\tau$', tau);

    figure;
    plot(1:ee, temp_corr(1,:), 1:ee, temp_corr(2,:), ...
        'Marker', '*', 'LineWidth', 2)
    set(gca, 'FontSize', 12)
    grid on
    title(titlestr1, 'Interpreter', 'latex', 'FontSize', 14)
    legend([names{pairs(ii,2)}, ' predicting ', names{pairs(ii,1)}], ...
        [names{pairs(ii,1)}, ' predicting ', names{pairs(ii,2)}])
    xlabel('Embedding dimension', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('Correlation ($\rho$)','Interpreter', 'latex', 'FontSize', 12)

    sprintf('Done with pair %d', ii)

    all_out.([names{pairs(ii,1)},'_',names{pairs(ii,2)}]) = temp_corr;
end