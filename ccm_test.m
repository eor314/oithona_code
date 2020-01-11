clear all
close all

%% Load data and get parameters
pt_oith = '/media/storage/learning_files/oithona_project/corrected_counts_091817.txt';

[raw, times] = load_oith(pt_oith, '03-11-2015 12:00:00 PM', '08-01-2015 12:00:00 AM');

% for ccm
tau = 24;
E = 5;

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

%test_oith = abs(hilbert(detrend((raw.oith - mean(raw.oith))./std(raw.oith))));
%test_egg = abs(hilbert(detrend((raw.egg - mean(raw.egg))./std(raw.egg))));

%% try the CCM
%[corrs , predY, predX, origY, origX] = sugi_CCM(filt.para, filt.oith, tau , E);
names = {'oith', 'egg', 'para'};
pairs = [1, 2; 1, 3; 2, 3];

all_out = struct();

for ii = 1:length(pairs)
    out = struct();
    for ee = 2:8
        temp_corr = zeros(2, 24);
        for tau = 1:24
            [corrs, ~] = sugi_CCM(filt.(names{pairs(ii,1)}), filt.(names{pairs(ii,2)}), ...
                tau, ee);
            temp_corr(:, tau) = corrs;
        end
        
        out.(sprintf('embed%d',ee)) = temp_corr;
        
        titlestr1 = sprintf([names{pairs(ii,2)}, ' predicted with ', ...
            names{pairs(ii,1)}, ' embedding %d'], ee);
        
        titlestr2 = sprintf([names{pairs(ii,1)}, ' predicted with ', ...
            names{pairs(ii,2)}, ' embedding %d'], ee);
        
        figure;
        subplot(2,1,1)
        plot(1:24, temp_corr(1,:))
        set(gca, 'FontSize', 12)
        grid on
        title(titlestr1, 'FontSize', 14)
        xlabel('$\tau$', 'Interpreter', 'latex', 'FontSize', 12)
        ylabel('Correlation', 'FontSize', 12)

        subplot(2,1,2)
        plot(1:24, temp_corr(2,:))
        set(gca, 'FontSize', 12)
        grid on
        title(titlestr2, 'FontSize', 14)
        xlabel('$\tau$', 'Interpreter', 'latex', 'FontSize', 12)
        ylabel('Correlation', 'FontSize', 12)

        sprintf('Done with embedding %d for pair %d', ee, ii)
    end
    all_out.([names{pairs(ii,1)},'_',names{pairs(ii,2)}]) = out;
end