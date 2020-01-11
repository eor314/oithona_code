%%%% Simple dynamic link model %%%%%

function [pred] = simple_dynamic_link(xx, deltat, tau)

%     simple dynamic link model for oithona project
%
%     :param xx: driver time-series for model (array)
%     :param deltat: desired timestep
%     :param tau: time lag (float)
%     :return pred: predicted time-series variable

pred(1) = 0;

for tt = 1:max(size(xx)) - 1
    pred(tt+1) = pred(tt) + (deltat/tau)*(xx(tt) - pred(tt));
end

end
