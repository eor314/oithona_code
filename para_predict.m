clear all 
close all

%%%%% Read data %%%%%
data = dlmread('/Users/Orenstein/Desktop/corrected_counts_091817.txt',',');

%%%%% Truncate zeros %%%%%
data = data(26:end-744, :);

%%%%% Make the corrected hourly counts and bin %%%%%
input1 = {'oith', data(1:end-1, 2) - data(1:end-1, 3) + data(1:end-1, 4)};
input2 = {'para', data(1:end-1, 8) - data(1:end-1, 9) + data(1:end-1, 10)};
input3 = {'egg', data(1:end-1, 5) - data(1:end-1, 6) + data(1:end-1, 7)};
time_frames = [1, 6, 24]; % time frames of interest

mm = 'run_mean'; % define method for averaging
%mm = 'block';

dd = bin_pack(time_frames, mm, input1, input2, input3);

%%
%%%%% loop over timeframes with simple model and plot %%%%%

names = fieldnames(dd);
%tau = 1:0.25:7; % define the range of time lags to compute
std_int = 24; % define the time interval of interest for the lags
%rr = zeros(max(size(names)), max(size(tau)));
rr = {};
pp = {};
pred_out = {};
titlestr = 'r^2 values, dynamic link at intervals: ';

%figure
%hold on
%grid on

comp = {'para', 'oith'};
titlestr2 = [comp{2}, ' driven by ', comp{1}, ' at intervals: '];

for ii = 1:max(size(names))
    
    % unpack the data structure
    %{
    oith_scaled = (dd.(names{ii}).oith - mean(dd.(names{ii}).oith)) / ...
        std(dd.(names{ii}).oith);
    para_scaled = (dd.(names{ii}).para - mean(dd.(names{ii}).para)) / ...
        std(dd.(names{ii}).para);
    %}
    obs_scaled = (dd.(names{ii}).(comp{1}) - mean(dd.(names{ii}).(comp{1}))) / ...
        std(dd.(names{ii}).(comp{1}));
    predict_scaled = (dd.(names{ii}).(comp{2}) - mean(dd.(names{ii}).(comp{2}))) / ...
        std(dd.(names{ii}).(comp{2}));
    % get the time step and compute the appropriate units for tau
    %ts = time_frames(ii); % time step between data points in hours
    %tau_unit = std_int/ts; % scale factor to make tau in units of std_int
    
    %tau_new = tau.*tau_unit; % recompute tau
    %tau_new = tau;
    
    %dt = ts/std_int; % set deltat 
    %dt = 1;
    
    dt = time_frames(ii)/std_int;
    %tau = (1/dt):(dt/4):(5/dt);
    tau = 1:(dt/4):(5/dt);
    r_ = zeros(size(tau));
    p_ = zeros(size(tau));
    
    % do the prediction
    for jj = 1:max(size(tau))
        preds = simple_dynamic_link(obs_scaled, dt, tau(jj));
        [r_temp, p_temp] = corrcoef(preds, predict_scaled); % get the r^2 value
        r_(jj) = r_temp(1,2);
        p_(jj) = p_temp(1,2);
        if jj > 1
            if r_(jj) > r_(jj-1)
                pred_best = preds;
            end
        else
            pred_best = preds;
        end
        clear temp
    end
    
    rr{ii} = r_;
    pp{ii} = p_;
    pred_out{ii} = pred_best;
    % generate plots
    
    % get datetime for axis label
    tt = 1:max(size(pred_best));
    startDate = datenum('03-11-2015 12:00:00 PM');
    endDate = datenum('08-01-2015 12:00:00 AM');
    xx = linspace(startDate, endDate, max(size(pred_best)));
    str = datestr(xx, 'mm-dd');
    
    % r^2 series
    figure
    titleout = [titlestr, ' ', names{ii}];
    plot(tau, r_, '-o', 'LineWidth', 2, 'MarkerSize', 8)
    grid on
    title(titleout, 'FontSize', 14)
    xlabel('time lag', 'FontSize', 12)
    ylabel('r^2', 'FontSize', 12)
    
    % predicted and actual
    figure
    titleout2 = [titlestr2, ' ', names{ii}];
    plot(xx, predict_scaled, xx, pred_best, 'LineWidth', 2)
    labs = get(gca,'XLim');
    set(gca, 'XTick', xx)
    axis tight
    datetick('x', 'mm-dd', 'keeplimits')
    grid on
    title(titleout2, 'FontSize', 14)
    legend('Measured', 'Predicted')
    xlabel('Time', 'FontSize', 12)
    ylabel('Number of ROIs', 'FontSize', 12)
end

%%
%%%%% loop over timeframes with more complex model and plot %%%%%
%{
names = fieldnames(dd);
%[tau1, tau2] = meshgrid(1:0.25:7); % define pairs of time lags to compute
rr = zeros(max(size(names)), size(tau1,1), size(tau1,2));
titlebase = 'r^2 values, dynamic link at intervals: ';
std_int = 24;

for ii = 1:max(size(names))
    
    % unpack the data structure
    oith_scaled = (dd.(names{ii}).oith - mean(dd.(names{ii}).oith)) / ...
        std(dd.(names{ii}).oith);
    para_scaled = (dd.(names{ii}).para - mean(dd.(names{ii}).para)) / ...
        std(dd.(names{ii}).para);
    
    dt = time_frames(ii)/std_int;
    [tau1, tau2] = meshgrid((1/dt):(dt/4):(5/dt));
    
    rr_temp = zeros(numel(tau1),1);
    % do the prediction
    for jj = 1:numel(tau1)
        preds = dynamic_link_posneg(oith_scaled, dt, tau1(jj), tau2(jj));
        temp = corrcoef(preds, para_scaled); % get the r^2 value
        rr_temp(jj) = temp(1, 2);
        clearvars temp
    end
    
    %reshape rr and save it out
    rr_out = reshape(rr_temp, size(tau1, 1), size(tau1, 2));
    rr(ii, :, :) = rr_out;
    
    % generate plots
    figure
    titlestr = [titlebase, ' ', names{ii}];
    surf(tau1, tau2, rr_out)
    title(titlestr, 'FontSize', 14)
    xlabel('tau1 (pos time lag)', 'FontSize', 12)
    ylabel('tau2 (neg time lag)', 'FontSize', 12)
    zlabel('r^2', 'FontSize', 12)

end
%}