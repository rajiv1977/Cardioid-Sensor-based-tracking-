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
function outIndex = deterministicR(inIndex,q)
% PURPOSE : Performs the resampling stage of the SIR
%           in order(number of samples) steps. It uses Kitagawa's
%           deterministic resampling algorithm.
% INPUTS  : - inIndex = Input particle indices.
%           - q = Normalised importance ratios.
% OUTPUTS : - outIndex = Resampled indices.

if nargin < 2, error('Not enough input arguments.'); end

[S,arb] = size(q);  % S = Number of particles.

% RESIDUAL RESAMPLING:
% ====================

N_babies= zeros(1,S);
u=zeros(1,S);

% generate the cumulative distribution
cumDist = cumsum(q');
aux=rand(1);
u=aux:1:(S-1+aux);
u=u./S;
j=1;
for i=1:S
   while (u(1,i)>cumDist(1,j))
      j=j+1;
   end
   N_babies(1,j)=N_babies(1,j)+1;
end

% COPY RESAMPLED TRAJECTORIES:  
% ============================
index=1;
for i=1:S
  if (N_babies(1,i)>0)
    for j=index:index+N_babies(1,i)-1
      outIndex(j) = inIndex(i);
    end;
  end;   
  index= index+N_babies(1,i);   
end












