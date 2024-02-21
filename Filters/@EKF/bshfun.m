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
function [y] = bshfun(x,u,t)
% PURPOSE : Measurement model function.
% INPUTS  : - x:  The evaluation point in the domain.
% OUTPUTS : - y:  The value of the function at x.

if nargin < 3, error('Not enough input arguments.'); end

[dim,np] = size(x);

y=zeros(2,np);

for j=1:np,

  r   = x(1,j);   % Risk free interest rate.
  sig = x(2,j);   % Volatility.
  S   = u(1);
  tm  = u(2);
  
  d1 = ( log(S) + (r+0.5*(sig^2))*tm ) / (sig * (tm^0.5));
  d2 = d1 - sig * (tm^0.5);
  % Compute call prices:
  y(1,j) = S*normcdf(d1) - exp(-r*tm)*normcdf(d2);
  % Compute put prices:
  y(2,j) = - S*normcdf(-d1) + exp(-r*tm)*normcdf(-d2);
 
end







