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
function [y] = ukf_bsffun(x,u,n,t)

% PURPOSE : Process model function for UKF.
% INPUTS  : - x:  The evaluation point in the domain.
%           - u:  control vector  
%           - n:  process noise
%           - t:  time index  
% OUTPUTS : - y:  The value of the function at x.

if nargin < 3, error('Not enough input arguments.'); end

y = bsffun(x,t) + n;









