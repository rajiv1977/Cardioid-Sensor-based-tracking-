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
% function fltr = EKF(sys, init_distr)
% EKF: class constructor for Extended Kalman Filter
%
% sys       : a SigObsSys object with Jacobians specified
% init_distr: initial distribution of the particles
% fltr      : the contructed filter
function fltr = EKF(sys, init_distr)

if nargin == 0
    fltr.f = [];
    fltr.h = [];
    fltr.Q = [];
    fltr.R = [];
    fltr.df = [];
    fltr.dh = [];
    fltr.x = [];
    fltr.P = [];
    fltr.phi = [];
    
    fltr = class(fltr, 'EKF');
elseif isa(sys, 'EKF')
    fltr = sys;
elseif nargin ~=2
    error('wrong number of arguments');
else
    if ~isa(sys, 'SigObsSys')
        error('wrong first argument: must be a SigObsSys object');
    end;
    [x, y, f, h, w_d, v_d, phi, df, dh] = get(sys);
    fltr.f = f;
    fltr.h = h;
    fltr.Q = cov(w_d);
    fltr.R = cov(v_d);
    fltr.df = df;
    fltr.dh = dh;
    fltr.x = mean(init_distr);
    fltr.P = cov(init_distr);
    fltr.phi = phi;
    
    fltr = class(fltr, 'EKF');
end;

end
