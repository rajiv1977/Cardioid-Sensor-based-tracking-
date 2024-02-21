%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Cardioid Sensor based tracking                                  %
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
% function obj =
% SigObsSys(fcn_sig,fcn_obs,distr_sig,distr_obs,x0,fcn_df, fcn_dh)
% SigObsSys: class constructor for a Signal-Observation system,
% which is a time-invariant dynamical system of the form:
%     x(k+1) = f(x(k)) + w(k)
%     y(k)   = h(x(k)) + v(k)
% fcn_sig is a function handle for f()
% fcn_obs is a function handle for h()
% distr_sig is a distribution object for w(k)
% distr_obs is a distribution object for v(k)
% x0 is an initial value for x(0)
% T0 is the sampling time
% fcn_df (optional) is Jacobian of f() with respect to x
% fcn_dh (optional) is Jacobian of h() with respect to x
function obj = SigObsSys(fcn_sig, fcn_obs, distr_sig, distr_obs, x0, T0, phi, fcn_df, fcn_dh)
if nargin == 0
    obj.f = [];
    obj.h = [];
    obj.w = [];
    obj.v = [];
    obj.x = [];
    obj.y = [];
    obj.df = [];
    obj.dh = [];
    obj.T = [];
    obj.phi = [];
    obj = class(obj,'SigObsSys');
elseif isa(fcn_sig, 'SigObsSys')
    obj = fcn_sig;
else
    switch nargin
        case 7,
            fcn_df = [];  fcn_dh = [];
        case 8,
            fcn_dh = [];
        case 9,
        otherwise,
            error('wrong number of arguments');
    end;
    
    obj.f = fcn_sig;
    obj.h = fcn_obs;
    obj.w = distr_sig;
    obj.v = distr_obs;
    obj.x = x0;
    obj.T = T0;
    obj.phi = phi;
    
    try
        % check for size in consistency
        len_x = length(obj.x);
        fx = obj.f(obj.x); % could fail, hence the try-catch
        len_f = length(fx);
        len_w = length(drawSamples(obj.w));
        obj.y = obj.h(obj.x, obj.T, obj.phi);  % could fail, hence the try-catch
        len_y = length(obj.y);
        len_v = length(drawSamples(obj.v));
        
        if len_x ~= len_w || len_x ~= len_f || len_y ~= len_v
            error('size inconsistency in SigObsSys constructor arguments');
        end;
        
    catch
        error('size inconsistency in SigObsSys constructor arguments');
    end
    
    obj.df = fcn_df;
    obj.dh = fcn_dh;
    obj = class(obj, 'SigObsSys');
end

end