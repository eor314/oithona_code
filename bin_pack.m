%%%%% Bin data and put into data structure %%%%%

function [dd] = bin_pack(timeframe, method, varargin)

%     Averages data over specified timeframes and returns a structure:
%
%     dd.{num_hours}_hour.{data_name} = [value1, value2, ... valueN]
%
%     Assumes original data is recorded in hourly intervals. 
%     
%     :param timeframes: array specifying the lengths of time frames (int)
%     :param method: indicates the desired method (string)
%     :param varargin: pairs of desired data name and associated array (cell)
%     :return dd: struct formated as dd.timeframe.dataname = binned array

% loop over timeframes
formatSpec = 'hour%d';
    
for tt = 1:max(size(timeframe))
    
    outstr = sprintf(formatSpec, timeframe(tt)); 
    
    % loop over pairs of data names and arrays
    for ii = 1:max(size(varargin))
        
        % read in pair
        in_cell = varargin{ii};
        
        % if the time frame is just 1, copy into struct
        if (timeframe(tt) == 1)
            dd.(outstr).(in_cell{1}) = in_cell{2};
        else
            if (strcmp(method, 'block'))
            % try to break into appropriate sized chunks
                try
                    out = mean(reshape(in_cell{2}, timeframe(tt), ...
                        max(size(in_cell{2}))/timeframe(tt)));

                    dd.(outstr).(in_cell{1}) = out;

                catch ME

                    if (strcmp(ME.identifier, 'MATLAB:getReshapeDims:notSameNumel'))
                        msg = ['Number of array elements cannot change when binning: ', ...
                            in_cell{1}, ' has ', num2str(max(size(in_cell{2}))), ...
                            ' elements and new reshaped array has dimsension', ...
                            num2str(tt), 'x', num2st(max(size(in_cell{2}))/tt), ...
                            '. Check input array size and desired bin width.'];
                        causeException = MException('MATLAB:bin_pack:dimensions',msg);
                        ME = addCause(ME, causeException);
                    end

                    rethrow(ME)
                end
                
            elseif (strcmp(method, 'run_mean'))
                % compute running mean in windows of size timeframe
                % (ignoring nans)
                dd.(outstr).(in_cell{1}) = movmean(in_cell{2}, timeframe(tt),...
                    'omitnan');
            else
                sprintf('Undefined averaging method')
            end
            
        end
    end
end
end