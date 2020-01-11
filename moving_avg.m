function out = moving_avg(win_size)

% make a moving average filter for data smoothing
% :param win_size: window size for filter
% :return out: 1-D filter

out = (1/win_size)*ones(1, win_size);
