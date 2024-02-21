%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Cardioid Sensor based tracking                                  %
%                            Auxiliary  Particle Filter                                     %
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
% input:
% fltr: the filter object
% y   : the current measurement
% output:
% fltr: the updated filter object

function fltr = update(fltr, y, ti)
N = length(fltr.w);
mu = fltr.p;
alpha = zeros(1, N);

for i = 1:N
    mu(:, i) = fltr.f(fltr.p(:, i));
    yp = fltr.h(mu(:, i), ti, fltr.phi);
    alpha(i) = density(fltr.v_d, y - yp);
end;

wa = fltr.w .* alpha;
sumwa = sum(wa);
if sumwa <= realmin
    error('adjustment multiplier is numerically zero in Auxiliary PF');
end;
wa = wa / sumwa;
outIndex = fcn_ResampSimp(wa);

w_smpl = drawSamples(fltr.w_d, N);

for i = 1:N
    j = outIndex(i);
    fltr.p(:, i) = mu(:, j) + w_smpl(:, i);
    y_pi = fltr.h(fltr.p(:, i), ti, fltr.phi);
    fltr.w(i) = density(fltr.v_d, y - y_pi) / alpha(j);
end;

sum_w = sum(fltr.w);
if sum_w <= realmin
    error('weights are numerically zero; change parameters or method.');
end;
fltr.w = fltr.w / sum_w;

end
