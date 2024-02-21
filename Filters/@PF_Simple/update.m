%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Cardioid Sensor based tracking                                  %
%                            Regular Particle Filter                                        %
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
function fltr = update(fltr, y, TT)
%fltr.r_fun
N = length(fltr.w);

% check whether resampling is chosen, and whether it's time to resample
if ~isempty(fltr.r_fun) && mod(fltr.counter, fltr.T) == 0
    if(N < fltr.thresh * fltr.N ) % too few (as a result of branch-kill), rejuvenate
        outIndex = fcn_ResampSys(fltr.w, fltr.N);
    else
        outIndex = fltr.r_fun(fltr.w);
    end;
    fltr.p = fltr.p(:, outIndex);
    N = length(outIndex);
    fltr.w = ones(1, N) / N;
end;

% internally keep track of when it's time to resample
fltr.counter = fltr.counter + 1;

w_smpl = drawSamples(fltr.w_d, N);
for i = 1:N
    
    % propagate particles
    fltr.p(:, i) = fltr.f(fltr.p(:, i)) + w_smpl(:, i);
    
    % noise-free y
    y_pi = fltr.h(fltr.p(:, i), TT, fltr.phi);
    
    % update weights
    fltr.w(i) = fltr.w(i) * density(fltr.v_d, y - y_pi);
    
end;

sum_w = sum(fltr.w);
if sum_w <= realmin
    error('weights are numerically zero; change parameters or method.');
end;

fltr.w = fltr.w / sum_w;

end
