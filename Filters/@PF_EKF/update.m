%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Cardioid Sensor based tracking                                  %
%                                     PF-EKF                                                %
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
function fltr = update(fltr, y, T)

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

for i = 1:N % for each particle
 
    x = fltr.p(:, i);
    xp = fltr.f(x);
    H = fltr.dh(xp, T, fltr.phi);
    yp = fltr.h(xp, T, fltr.phi);
    P = inv(fltr.invQ + H' * fltr.invR * H);
    P = (P + P') / 2; % P can become numerically asymmetric
       
    mu = P * (fltr.invQ * xp + H * fltr.invR * (y - yp + H' * xp));
    
    gauss = GaussianDistr(mu, P); % proposal distribution
    
    % propagate particles
    fltr.p(:, i) = drawSamples(gauss, 1);
    % mean of y
    y_pi = fltr.h(fltr.p(:, i), T, fltr.phi);
    % from measurement model
    likelihood = density(fltr.v_d, y - y_pi);
    % from process model
    prior = density(fltr.w_d, (fltr.p(:,i) - xp));
 
    % from the proposal
    proposal = density(gauss, fltr.p(:,i));

    fltr.w(i) = fltr.w(i) * likelihood * prior / proposal;
    
end;
  
sum_w = sum(fltr.w);
if sum_w <= realmin^2
    error('weights are numerically zero; change parameters or method.');
end;

fltr.w = fltr.w / sum_w;

end
