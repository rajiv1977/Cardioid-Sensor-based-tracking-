%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Cardioid Sensor based tracking                                  %
%                                        EKF                                                %
%                     Copyright @2015_DRDC, version 01_02112015                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv,  and B.Balaji                                      %
%          Defence R&D Canada, 3701 Carling Avenue, Ottawa, ON, K1A 0Z4, Canada.            %
%             rajiv.sithiravel@gmail.com and Bhashyam.Balaji@drdc-rddc.gc.ca                %
%                                                                                           %
%                                   T.Kirubarajan                                           %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                                 kiruba@mcmaster.ca                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fltr = update(fltr, y)
%
% INPUT:
%        fltr: the filter
%           y: the measurement
% OUTPUT:
%        fltr: updated filter

function fltr = update(fltr, y, T)

try
    A = fltr.df(fltr.x);
    P = A * fltr.P * A' + fltr.Q;
    x = fltr.f(fltr.x);
    H = fltr.dh(x, T, fltr.phi);
    K = P * H' * inv(H * P * H' + fltr.R);
    yhat = fltr.h(x, T, fltr.phi);
    fltr.x = x + K * (y - yhat);
    fltr.P = P - K * H * P;
catch
    error('EKF initialized incorrectly with a SigObsSys object');
end;
end

