%%%%% Two tau dynamic link model %%%%%

function [pred] = dynamic_link_posneg(xx, deltat, tau1, tau2)

%     dynamic link model that applies a different time lag depending on 
%     direction of change of driver variable
%
%     :param xx: driver time-series for model (array)
%     :param deltat: desired timestep
%     :param tau1: time lag if driver is positive change (float)
%     :param tau2: time lag if driver is negative change (float)
%     :return pred: predicted time-series variable

pred(1) = 0;
for tt = 1:max(size(xx)) - 1
    
    % assume that change is positive in first time point of driver
    if (tt == 1)
        pred(tt+1) = pred(tt) + (deltat/tau1)*(xx(tt) - pred(tt));
    else
        if (sign(xx(tt) - xx(tt-1)) == 1)
            pred(tt+1) = pred(tt) + (deltat/tau1)*(xx(tt) - pred(tt));
        else
            pred(tt+1) = pred(tt) + (deltat/tau2)*(xx(tt) - pred(tt));
        end
    end
end
end