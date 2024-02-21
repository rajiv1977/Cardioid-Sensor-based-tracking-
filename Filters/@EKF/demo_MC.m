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
% PURPOSE : Demonstrate the differences between the following filters on the same problem:
%           
%           1) Extended Kalman Filter  (EKF)
%           2) Unscented Kalman Filter (UKF)
%           3) Particle Filter         (PF)
%           4) PF with EKF proposal    (PFEKF)
%           5) PF with UKF proposal    (PFUKF)

clear all;
clc;
echo off;
path('./ukf',path);

% INITIALISATION AND PARAMETERS:
% ==============================

no_of_runs = 100;            % number of experiments to generate statistical
                            % averages
doPlot = 0;                 % 1 plot online. 0 = only plot at the end.
sigma =  1e-5;              % Variance of the Gaussian measurement noise.
g1 = 3;                     % Paramater of Gamma transition prior.
g2 = 2;                     % Parameter of Gamman transition prior.
                            % Thus mean = 3/2 and var = 3/4.
T = 60;                     % Number of time steps.
R = 1e-5;                   % EKF's measurement noise variance. 
Q = 3/4;                    % EKF's process noise variance.
P0 = 3/4;                   % EKF's initial variance of the states.

N = 200;                     % Number of particles.
resamplingScheme = 1;       % The possible choices are
                            % systematic sampling (2),
                            % residual (1)
                            % and multinomial (3). 
                            % They're all O(N) algorithms. 

Q_pfekf = 10*3/4;
R_pfekf = 1e-1;

Q_pfukf = 2*3/4;
R_pfukf = 1e-1;
			    
alpha = 1;                  % UKF : point scaling parameter
beta  = 0;                  % UKF : scaling parameter for higher order terms of Taylor series expansion 
kappa = 2;                  % UKF : sigma point selection scaling parameter (best to leave this = 0)

%**************************************************************************************

% SETUP BUFFERS TO STORE PERFORMANCE RESULTS
% ==========================================

rmsError_ekf      = zeros(1,no_of_runs);
rmsError_ukf      = zeros(1,no_of_runs);
rmsError_pf       = zeros(1,no_of_runs);
rmsError_pfMC     = zeros(1,no_of_runs);
rmsError_pfekf    = zeros(1,no_of_runs);
rmsError_pfekfMC  = zeros(1,no_of_runs);
rmsError_pfukf    = zeros(1,no_of_runs);
rmsError_pfukfMC  = zeros(1,no_of_runs);

time_pf       = zeros(1,no_of_runs);     
time_pfMC     = zeros(1,no_of_runs);
time_pfekf    = zeros(1,no_of_runs);
time_pfekfMC  = zeros(1,no_of_runs);
time_pfukf    = zeros(1,no_of_runs);
time_pfukfMC  = zeros(1,no_of_runs);

%**************************************************************************************

% MAIN LOOP

for j=1:no_of_runs,

  rand('state',sum(100*clock));   % Shuffle the pack!
  randn('state',sum(100*clock));   % Shuffle the pack!
  

% GENERATE THE DATA:
% ==================

x = zeros(T,1);
y = zeros(T,1);
processNoise = zeros(T,1);
measureNoise = zeros(T,1);
x(1) = 1;                         % Initial state.
for t=2:T
  processNoise(t) = gengamma(g1,g2);  
  measureNoise(t) = sqrt(sigma)*randn(1,1);    
  x(t) = feval('ffun',x(t-1),t) +processNoise(t);     % Gamma transition prior.  
  y(t) = feval('hfun',x(t),t) + measureNoise(t);      % Gaussian likelihood.
end;  

% PLOT THE GENERATED DATA:
% ========================
figure(1)
clf;
plot(1:T,x,'r',1:T,y,'b');
ylabel('Data','fontsize',15);
xlabel('Time','fontsize',15);
legend('States (x)','Observations(y)');

%%%%%%%%%%%%%%%  PERFORM EKF and UKF ESTIMATION  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  ==============================  %%%%%%%%%%%%%%%%%%%%%

% INITIALISATION:
% ==============
mu_ekf = ones(T,1);     % EKF estimate of the mean of the states.
P_ekf = P0*ones(T,1);   % EKF estimate of the variance of the states.

mu_ukf = mu_ekf;        % UKF estimate of the mean of the states.
P_ukf = P_ekf;          % UKF estimate of the variance of the states.

yPred = ones(T,1);      % One-step-ahead predicted values of y.
mu_ekfPred = ones(T,1); % EKF O-s-a estimate of the mean of the states.
PPred = ones(T,1);      % EKF O-s-a estimate of the variance of the states.

disp(' ');

for t=2:T,    
  fprintf('run = %i / %i :  EKF & UKF : t = %i / %i  \r',j,no_of_runs,t,T);
  fprintf('\n')
  
  % PREDICTION STEP:
  % ================ 
  mu_ekfPred(t) = feval('ffun',mu_ekf(t-1),t);
  Jx = 0.5;                             % Jacobian for ffun.
  PPred(t) = Q + Jx*P_ekf(t-1)*Jx'; 
  
  % CORRECTION STEP:
  % ================
  yPred(t) = feval('hfun',mu_ekfPred(t),t);
  if t<=30,
    Jy = 2*0.2*mu_ekfPred(t);                 % Jacobian for hfun.
  else
    Jy = 0.5;
  %  Jy = cos(mu_ekfPred(t))/2;
  %   Jy = 2*mu_ekfPred(t)/4;                 % Jacobian for hfun. 
  end;
  M = R + Jy*PPred(t)*Jy';                 % Innovations covariance.
  K = PPred(t)*Jy'*inv(M);                 % Kalman gain.
  mu_ekf(t) = mu_ekfPred(t) + K*(y(t)-yPred(t));
  P_ekf(t) = PPred(t) - K*Jy*PPred(t);
  
  % Full Unscented Kalman Filter step
  % =================================
  [mu_ukf(t),P_ukf(t)]=ukf(mu_ukf(t-1),P_ukf(t-1),[],Q,'ukf_ffun',y(t),R,'ukf_hfun',t,alpha,beta,kappa);  
  
end;   % End of t loop.



%%%%%%%%%%%%%%%  PERFORM SEQUENTIAL MONTE CARLO  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  ==============================  %%%%%%%%%%%%%%%%%%%%%

% INITIALISATION:
% ==============
xparticle_pf = ones(T,N);        % These are the particles for the estimate
                                 % of x. Note that there's no need to store
                                 % them for all t. We're only doing this to
                                 % show you all the nice plots at the end.
xparticlePred_pf = ones(T,N);    % One-step-ahead predicted values of the states.
yPred_pf = ones(T,N);            % One-step-ahead predicted values of y.
w = ones(T,N);                   % Importance weights.

disp(' ');
 
tic;                             % Initialize timer for benchmarking

for t=2:T,    
  fprintf('run = %i / %i :  PF : t = %i / %i  \r',j,no_of_runs,t,T);
  fprintf('\n')
  
  % PREDICTION STEP:
  % ================ 
  % We use the transition prior as proposal.
  for i=1:N,
    xparticlePred_pf(t,i) = feval('ffun',xparticle_pf(t-1,i),t) + gengamma(g1,g2);   
  end;

  % EVALUATE IMPORTANCE WEIGHTS:
  % ============================
  % For our choice of proposal, the importance weights are give by:  
  for i=1:N,
    yPred_pf(t,i) = feval('hfun',xparticlePred_pf(t,i),t);        
    lik = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-yPred_pf(t,i))^(2))) ...
	  + 1e-99; % Deal with ill-conditioning.
    w(t,i) = lik;    
  end;  
  w(t,:) = w(t,:)./sum(w(t,:));                % Normalise the weights.
  
  % SELECTION STEP:
  % ===============
  % Here, we give you the choice to try three different types of
  % resampling algorithms. Note that the code for these algorithms
  % applies to any problem!
  if resamplingScheme == 1
    outIndex = residualR(1:N,w(t,:)');        % Residual resampling.
  elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w(t,:)');      % Systematic resampling.
  else  
    outIndex = multinomialR(1:N,w(t,:)');     % Multinomial resampling.  
  end;
  xparticle_pf(t,:) = xparticlePred_pf(t,outIndex); % Keep particles with
                                                    % resampled indices.
end;   % End of t loop.

time_pf(j) = toc;    % How long did this take?


%%%%%%%%%%%%%%  PERFORM SEQUENTIAL MONTE CARLO WITH MCMC  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  ========================================  %%%%%%%%%%%%%%%%

% INITIALISATION:
% ==============
xparticle_pfMC = ones(T,N);      % These are the particles for the estimate
                                 % of x. Note that there's no need to store
                                 % them for all t. We're only doing this to
                                 % show you all the nice plots at the end.
xparticlePred_pfMC = ones(T,N);  % One-step-ahead predicted values of the states.
yPred_pfMC = ones(T,N);          % One-step-ahead predicted values of y.
w = ones(T,N);                   % Importance weights.
previousXMC = ones(T,N);         % Particles at the previous time step. 
previousXResMC = ones(T,N);      % Resampled previousX.

disp(' ');
 
tic;                             % Initialize timer for benchmarking

for t=2:T,    
  fprintf('run = %i / %i :  PF-MCMC : t = %i / %i  \r',j,no_of_runs,t,T);
  fprintf('\n')
  
  % PREDICTION STEP:
  % ================ 
  % We use the transition prior as proposal.
  for i=1:N,
    xparticlePred_pfMC(t,i) = feval('ffun',xparticle_pfMC(t-1,i),t) + gengamma(g1,g2);   
  end;
  previousXMC(t,:) = xparticle_pfMC(t-1,:);  % Store the particles at t-1. 


  % EVALUATE IMPORTANCE WEIGHTS:
  % ============================
  % For our choice of proposal, the importance weights are give by:  
  for i=1:N,
    yPred_pfMC(t,i) = feval('hfun',xparticlePred_pfMC(t,i),t);        
    lik = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-yPred_pfMC(t,i))^(2))) ...
	  + 1e-99; % Deal with ill-conditioning.
    w(t,i) = lik;    
  end;  
  w(t,:) = w(t,:)./sum(w(t,:));                % Normalise the weights.
  
  % SELECTION STEP:
  % ===============
  % Here, we give you the choice to try three different types of
  % resampling algorithms. Note that the code for these algorithms
  % applies to any problem!
  if resamplingScheme == 1
    outIndex = residualR(1:N,w(t,:)');        % Residual resampling.
  elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w(t,:)');      % Systematic resampling.
  else  
    outIndex = multinomialR(1:N,w(t,:)');     % Multinomial resampling.  
  end;
  xparticle_pfMC(t,:) = xparticlePred_pfMC(t,outIndex); % Keep particles with
                                                        % resampled
                                                        % indices.
  previousXResMC(t,:) = previousXMC(t,outIndex);  % Resample particles
                                                  % at t-1.
  
  % METROPOLIS-HASTINGS STEP:
  % ========================
  u=rand(N,1); 
  accepted=0;
  rejected=0;
  for i=1:N,   
    xProp = feval('ffun',previousXResMC(t,i),t) + gengamma(g1,g2);   
    mProp = feval('hfun',xProp,t);        
    likProp = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-mProp)^(2))) + 1e-99;     
    m = feval('hfun',xparticle_pfMC(t,i),t);        
    lik = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-m)^(2))) + 1e-99;     
    acceptance = min(1,likProp/lik);
    if u(i,1) <= acceptance 
      xparticle_pfMC(t,i) = xProp;
      accepted=accepted+1;
    else
      xparticle_pfMC(t,i) = xparticle_pfMC(t,i); 
      rejected=rejected+1;
    end;
  end;
  
  
end;   % End of t loop.

time_pfMC(j) = toc;    % How long did this take?


%%%%%%%%%%%%%%%  PERFORM SEQUENTIAL MONTE CARLO  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  ======== EKF proposal ========  %%%%%%%%%%%%%%%%%%%%%

% INITIALISATION:
% ==============
xparticle_pfekf = ones(T,N);        % These are the particles for the estimate
                                    % of x. Note that there's no need to store
                                    % them for all t. We're only doing this to
                                    % show you all the nice plots at the end.
Pparticle_pfekf = P0*ones(T,N);     % Particles for the covariance of x.
xparticlePred_pfekf = ones(T,N);    % One-step-ahead predicted values of the states.
PparticlePred_pfekf = ones(T,N);    % One-step-ahead predicted values of P.
yPred_pfekf = ones(T,N);            % One-step-ahead predicted values of y.
w = ones(T,N);                      % Importance weights.
muPred_pfekf = ones(T,1);           % EKF O-s-a estimate of the mean of the states.
PPred_pfekf = ones(T,1);            % EKF O-s-a estimate of the variance of the states.
mu_pfekf = ones(T,1);               % EKF estimate of the mean of the states.
P_pfekf = P0*ones(T,1);             % EKF estimate of the variance of the states.

disp(' ');

tic;                                % Initialize timer for benchmarking

for t=2:T,    
  fprintf('run = %i / %i :  PF-EKF : t = %i / %i  \r',j,no_of_runs,t,T);
  fprintf('\n')
  
  % PREDICTION STEP:
  % ================ 
  % We use the EKF as proposal.
  for i=1:N,
    muPred_pfekf(t) = feval('ffun',xparticle_pfekf(t-1,i),t);
    Jx = 0.5;                                 % Jacobian for ffun.
    PPred_pfekf(t) = Q_pfekf + Jx*Pparticle_pfekf(t-1,i)*Jx'; 
    yPredTmp = feval('hfun',muPred_pfekf(t),t);
    if t<=30,
      Jy = 2*0.2*muPred_pfekf(t);                     % Jacobian for hfun.
    else
      Jy = 0.5;
    end;
    M = R_pfekf + Jy*PPred_pfekf(t)*Jy';                  % Innovations covariance.
    K = PPred_pfekf(t)*Jy'*inv(M);                  % Kalman gain.
    mu_pfekf(t,i) = muPred_pfekf(t) + K*(y(t)-yPredTmp); % Mean of proposal.
    P_pfekf(t) = PPred_pfekf(t) - K*Jy*PPred_pfekf(t);          % Variance of proposal.
    xparticlePred_pfekf(t,i) = mu_pfekf(t,i) + sqrtm(P_pfekf(t))*randn(1,1);
    PparticlePred_pfekf(t,i) = P_pfekf(t);
  end;

  % EVALUATE IMPORTANCE WEIGHTS:
  % ============================
  % For our choice of proposal, the importance weights are give by:  
  for i=1:N,
    yPred_pfekf(t,i) = feval('hfun',xparticlePred_pfekf(t,i),t);        
    lik = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-yPred_pfekf(t,i))^(2)))+1e-99;
    prior = ((xparticlePred_pfekf(t,i)-xparticle_pfekf(t-1,i))^(g1-1)) ...
		 * exp(-g2*(xparticlePred_pfekf(t,i)-xparticle_pfekf(t-1,i)));
    proposal = inv(sqrt(PparticlePred_pfekf(t,i))) * ...
	       exp(-0.5*inv(PparticlePred_pfekf(t,i)) *((xparticlePred_pfekf(t,i)-mu_pfekf(t,i))^(2)));
    w(t,i) = lik*prior/proposal;      
  end;  
  w(t,:) = w(t,:)./sum(w(t,:));                % Normalise the weights.
  
  % SELECTION STEP:
  % ===============
  % Here, we give you the choice to try three different types of
  % resampling algorithms. Note that the code for these algorithms
  % applies to any problem!
  if resamplingScheme == 1
    outIndex = residualR(1:N,w(t,:)');        % Residual resampling.
  elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w(t,:)');      % Systematic resampling.
  else  
    outIndex = multinomialR(1:N,w(t,:)');     % Multinomial resampling.  
  end;
  xparticle_pfekf(t,:) = xparticlePred_pfekf(t,outIndex); % Keep particles with
                                              % resampled indices.
  Pparticle_pfekf(t,:) = PparticlePred_pfekf(t,outIndex);  
  
end;   % End of t loop.

time_pfekf(j) = toc;

%%%%%%%%%%%%%%  PERFORM SEQUENTIAL MONTE CARLO WITH MCMC  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  ======== EKF proposal ==================  %%%%%%%%%%%%%%%%

% INITIALISATION:
% ==============
xparticle_pfekfMC = ones(T,N);        % These are the particles for the estimate
                                      % of x. Note that there's no need to store
                                      % them for all t. We're only doing this to
                                      % show you all the nice plots at the end.
Pparticle_pfekfMC = P0*ones(T,N);     % Particles for the covariance of x.
xparticlePred_pfekfMC = ones(T,N);    % One-step-ahead predicted values of the states.
PparticlePred_pfekfMC = ones(T,N);    % One-step-ahead predicted values of P.
yPred_pfekfMC = ones(T,N);            % One-step-ahead predicted values of y.
w = ones(T,N);                        % Importance weights.
muPred_pfekfMC = ones(T,1);           % EKF O-s-a estimate of the mean of the states.
PPred_pfekfMC = ones(T,1);            % EKF O-s-a estimate of the variance of the states.
mu_pfekfMC = ones(T,1);               % EKF estimate of the mean of the states.
P_pfekfMC = P0*ones(T,1);             % EKF estimate of the variance of the states.

previousXekfMC = ones(T,N);           % Particles at the previous time step. 
previousXResekfMC = ones(T,N);        % Resampled previousX.
previousPekfMC = ones(T,N);           % Covariance particles at the previous time step. 
previousPResekfMC = ones(T,N);        % Resampled previousP.


disp(' ');

tic;                                % Initialize timer for benchmarking

for t=2:T,    
  fprintf('run = %i / %i :  PF-EKF-MCMC : t = %i / %i  \r',j,no_of_runs,t,T);
  fprintf('\n')
  
  % PREDICTION STEP:
  % ================ 
  % We use the EKF as proposal.
  for i=1:N,
    muPred_pfekfMC(t) = feval('ffun',xparticle_pfekfMC(t-1,i),t);
    Jx = 0.5;                                 % Jacobian for ffun.
    PPred_pfekfMC(t) = Q_pfekf + Jx*Pparticle_pfekfMC(t-1,i)*Jx'; 
    yPredTmp = feval('hfun',muPred_pfekfMC(t),t);
    if t<=30,
      Jy = 2*0.2*muPred_pfekfMC(t);                     % Jacobian for hfun.
    else
      Jy = 0.5;
    end;
    M = R_pfekf + Jy*PPred_pfekfMC(t)*Jy';                  % Innovations covariance.
    K = PPred_pfekfMC(t)*Jy'*inv(M);                  % Kalman gain.
    mu_pfekfMC(t,i) = muPred_pfekfMC(t) + K*(y(t)-yPredTmp); % Mean of proposal.
    P_pfekfMC(t) = PPred_pfekfMC(t) - K*Jy*PPred_pfekfMC(t);          % Variance of proposal.
    xparticlePred_pfekfMC(t,i) = mu_pfekfMC(t,i) + sqrtm(P_pfekfMC(t))*randn(1,1);
    PparticlePred_pfekfMC(t,i) = P_pfekfMC(t);
  end;

  previousXekfMC(t,:) = xparticle_pfekfMC(t-1,:);  % Store the particles at t-1. 
  previousPekfMC(t,:) = Pparticle_pfekfMC(t-1,:);  % Store the particles at t-1. 
  
  
  % EVALUATE IMPORTANCE WEIGHTS:
  % ============================
  % For our choice of proposal, the importance weights are give by:  
  for i=1:N,
    yPred_pfekfMC(t,i) = feval('hfun',xparticlePred_pfekfMC(t,i),t);        
    lik = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-yPred_pfekfMC(t,i))^(2)))+1e-99;
    prior = ((xparticlePred_pfekfMC(t,i)-xparticle_pfekfMC(t-1,i))^(g1-1)) ...
		 * exp(-g2*(xparticlePred_pfekfMC(t,i)-xparticle_pfekfMC(t-1,i)));
    proposal = inv(sqrt(PparticlePred_pfekfMC(t,i))) * ...
	       exp(-0.5*inv(PparticlePred_pfekfMC(t,i)) *((xparticlePred_pfekfMC(t,i)-mu_pfekfMC(t,i))^(2)));
    w(t,i) = lik*prior/proposal;      
  end;  
  w(t,:) = w(t,:)./sum(w(t,:));                % Normalise the weights.
  
  % SELECTION STEP:
  % ===============
  % Here, we give you the choice to try three different types of
  % resampling algorithms. Note that the code for these algorithms
  % applies to any problem!
  if resamplingScheme == 1
    outIndex = residualR(1:N,w(t,:)');        % Residual resampling.
  elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w(t,:)');      % Systematic resampling.
  else  
    outIndex = multinomialR(1:N,w(t,:)');     % Multinomial resampling.  
  end;
  xparticle_pfekfMC(t,:) = xparticlePred_pfekfMC(t,outIndex); % Keep particles with
                                                         % resampled indices.
  Pparticle_pfekfMC(t,:) = PparticlePred_pfekfMC(t,outIndex);  
  previousXResekfMC(t,:) = previousXekfMC(t,outIndex);  % Resample particles
                                                        % at t-1.
  previousPResekfMC(t,:) = previousPekfMC(t,outIndex);  % Resample particles
                                                        % at t-1.
   
  % METROPOLIS-HASTINGS STEP:
  % ========================
  u=rand(N,1); 
  accepted=0;
  rejected=0;
  for i=1:N,   
    muPred_ekfMCMC = feval('ffun',previousXResekfMC(t,i),t);
    Jx = 0.5;                                     % Jacobian for ffun.
    PPred_ekfMCMC = Q_pfekf + Jx*previousPResekfMC(t,i)*Jx'; 
    yPredTmp = feval('hfun',muPred_ekfMCMC,t);
    if t<=30,
      Jy = 2*0.2*muPred_ekfMCMC;                     % Jacobian for hfun.
    else
      Jy = 0.5;
    end;
    M = R_pfekf + Jy*PPred_ekfMCMC*Jy';                  % Innovations covariance.
    K = PPred_ekfMCMC*Jy'*inv(M);                  % Kalman gain.
    muProp = muPred_ekfMCMC + K*(y(t)-yPredTmp);   % Mean of proposal.
    PProp = PPred_ekfMCMC - K*Jy*PPred_ekfMCMC;          % Variance of proposal.
    xparticleProp = muProp + sqrtm(PProp)*randn(1,1);
    PparticleProp = PProp;   
    
    mProp = feval('hfun',xparticleProp,t);        
    likProp = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-mProp)^(2)))+1e-99;
    priorProp = ((xparticleProp-previousXResekfMC(t,i))^(g1-1)) ...
		 * exp(-g2*(xparticleProp-previousXResekfMC(t,i)));
    proposalProp = inv(sqrt(PparticleProp)) * ...
	       exp(-0.5*inv(PparticleProp) *( ...
					      (xparticleProp-muProp)^(2)));
    m = feval('hfun',xparticle_pfekfMC(t,i),t);        
    lik = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-m)^(2)))+1e-99;
    prior = ((xparticle_pfekfMC(t,i)-previousXResekfMC(t,i))^(g1-1)) ...
		 * exp(-g2*(xparticle_pfekfMC(t,i)-previousXResekfMC(t,i)));
    proposal = inv(sqrt(Pparticle_pfekfMC(t,i))) * ...
	       exp(-0.5*inv(Pparticle_pfekfMC(t,i)) *((xparticle_pfekfMC(t,i)-muProp)^(2)));
    ratio = (likProp*priorProp*proposal)/(lik*prior*proposalProp);
    acceptance = min(1,ratio);
    if u(i,1) <= acceptance 
      xparticle_pfekfMC(t,i) = xparticleProp;
      Pparticle_pfekfMC(t,i) = PparticleProp;
      accepted=accepted+1;
    else
      xparticle_pfekfMC(t,i) = xparticle_pfekfMC(t,i); 
      Pparticle_pfekfMC(t,i) = Pparticle_pfekfMC(t,i);  
      rejected=rejected+1;
    end;
  end;   % End of MCMC loop.
end;   % End of t loop.

time_pfekfMC(j) = toc;


%%%%%%%%%%%%%%%  PERFORM SEQUENTIAL MONTE CARLO  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  ======== UKF proposal ========  %%%%%%%%%%%%%%%%%%%%%

% INITIALISATION:
% ==============
xparticle_pfukf = ones(T,N);        % These are the particles for the estimate
                                    % of x. Note that there's no need to store
                                    % them for all t. We're only doing this to
                                    % show you all the nice plots at the end.
Pparticle_pfukf = P0*ones(T,N);     % Particles for the covariance of x.
xparticlePred_pfukf = ones(T,N);    % One-step-ahead predicted values of the states.
PparticlePred_pfukf = ones(T,N);    % One-step-ahead predicted values of P.
yPred_pfukf = ones(T,N);            % One-step-ahead predicted values of y.
w = ones(T,N);                      % Importance weights.
mu_pfukf = ones(T,1);               % EKF estimate of the mean of the states.

error=0;

disp(' ');

tic;

for t=2:T,    
  fprintf('run = %i / %i :  PF-UKF : t = %i / %i  \r',j,no_of_runs,t,T);
  fprintf('\n')
  
  % PREDICTION STEP:
  % ================ 
  % We use the UKF as proposal.
  for i=1:N,
    % Call Unscented Kalman Filter
    [mu_pfukf(t,i),PparticlePred_pfukf(t,i)]=ukf(xparticle_pfukf(t-1,i),Pparticle_pfukf(t-1,i),[],Q_pfukf,'ukf_ffun',y(t),R_pfukf,'ukf_hfun',t,alpha,beta,kappa);
    xparticlePred_pfukf(t,i) = mu_pfukf(t,i) + sqrtm(PparticlePred_pfukf(t,i))*randn(1,1);
  end;

  % EVALUATE IMPORTANCE WEIGHTS:
  % ============================
  % For our choice of proposal, the importance weights are give by:  
  for i=1:N,
    yPred_pfukf(t,i) = feval('hfun',xparticlePred_pfukf(t,i),t);        
    lik = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-yPred_pfukf(t,i))^(2)))+1e-99;
    prior = ((xparticlePred_pfukf(t,i)-xparticle_pfukf(t-1,i))^(g1-1)) ...
		 * exp(-g2*(xparticlePred_pfukf(t,i)-xparticle_pfukf(t-1,i)));
    proposal = inv(sqrt(PparticlePred_pfukf(t,i))) * ...
	       exp(-0.5*inv(PparticlePred_pfukf(t,i)) *((xparticlePred_pfukf(t,i)-mu_pfukf(t,i))^(2)));
    w(t,i) = lik*prior/proposal;      
  end;  
  w(t,:) = w(t,:)./sum(w(t,:));                % Normalise the weights.
  
  % SELECTION STEP:
  % ===============
  % Here, we give you the choice to try three different types of
  % resampling algorithms. Note that the code for these algorithms
  % applies to any problem!
  if resamplingScheme == 1
    outIndex = residualR(1:N,w(t,:)');        % Residual resampling.
  elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w(t,:)');      % Systematic resampling.
  else  
    outIndex = multinomialR(1:N,w(t,:)');     % Multinomial resampling.  
  end;
  xparticle_pfukf(t,:) = xparticlePred_pfukf(t,outIndex); % Keep particles with
                                              % resampled indices.
  Pparticle_pfukf(t,:) = PparticlePred_pfukf(t,outIndex);  
  
end;   % End of t loop.

time_pfukf(j) = toc;





%%%%%%%%%%%%%%  PERFORM SEQUENTIAL MONTE CARLO WITH MCMC  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  ============= UKF proposal =============  %%%%%%%%%%%%%%%%

% INITIALISATION:
% ==============
xparticle_pfukfMC = ones(T,N);        % These are the particles for the estimate
                                      % of x. Note that there's no need to store
                                      % them for all t. We're only doing this to
                                      % show you all the nice plots at the end.
Pparticle_pfukfMC = P0*ones(T,N);     % Particles for the covariance of x.
xparticlePred_pfukfMC = ones(T,N);    % One-step-ahead predicted values of the states.
PparticlePred_pfukfMC = ones(T,N);    % One-step-ahead predicted values of P.
yPred_pfukfMC = ones(T,N);            % One-step-ahead predicted values of y.
w = ones(T,N);                        % Importance weights.
mu_pfukfMC = ones(T,1);               % EKF estimate of the mean of the states.

previousXukfMC = ones(T,N);           % Particles at the previous time step. 
previousXResukfMC = ones(T,N);        % Resampled previousX.
previousPukfMC = ones(T,N);           % Covariance particles at the previous time step. 
previousPResukfMC = ones(T,N);        % Resampled previousP.

error=0;

disp(' ');

tic;

for t=2:T,    
  fprintf('run = %i / %i :  PF-UKF-MCMC : t = %i / %i  \r',j,no_of_runs,t,T);
  fprintf('\n')
  
  % PREDICTION STEP:
  % ================ 
  % We use the UKF as proposal.
  for i=1:N,
    % Call Unscented Kalman Filter
    [mu_pfukfMC(t,i),PparticlePred_pfukfMC(t,i)]=ukf(xparticle_pfukfMC(t-1,i),Pparticle_pfukfMC(t-1,i),[],Q_pfukf,'ukf_ffun',y(t),R_pfukf,'ukf_hfun',t,alpha,beta,kappa);
    xparticlePred_pfukfMC(t,i) = mu_pfukfMC(t,i) + sqrtm(PparticlePred_pfukfMC(t,i))*randn(1,1);
  end;
  
  previousXukfMC(t,:) = xparticle_pfukfMC(t-1,:);  % Store the particles at t-1. 
  previousPukfMC(t,:) = Pparticle_pfukfMC(t-1,:);  % Store the particles at t-1. 
  
   

  % EVALUATE IMPORTANCE WEIGHTS:
  % ============================
  % For our choice of proposal, the importance weights are give by:  
  for i=1:N,
    yPred_pfukfMC(t,i) = feval('hfun',xparticlePred_pfukfMC(t,i),t);        
    lik = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-yPred_pfukfMC(t,i))^(2)))+1e-99;
    prior = ((xparticlePred_pfukfMC(t,i)-xparticle_pfukfMC(t-1,i))^(g1-1)) ...
		 * exp(-g2*(xparticlePred_pfukfMC(t,i)-xparticle_pfukfMC(t-1,i)));
    proposal = inv(sqrt(PparticlePred_pfukfMC(t,i))) * ...
	       exp(-0.5*inv(PparticlePred_pfukfMC(t,i)) *((xparticlePred_pfukfMC(t,i)-mu_pfukfMC(t,i))^(2)));
    w(t,i) = lik*prior/proposal;      
  end;  
  w(t,:) = w(t,:)./sum(w(t,:));                % Normalise the weights.
  
  % SELECTION STEP:
  % ===============
  % Here, we give you the choice to try three different types of
  % resampling algorithms. Note that the code for these algorithms
  % applies to any problem!
  if resamplingScheme == 1
    outIndex = residualR(1:N,w(t,:)');        % Residual resampling.
  elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w(t,:)');      % Systematic resampling.
  else  
    outIndex = multinomialR(1:N,w(t,:)');     % Multinomial resampling.  
  end;
  xparticle_pfukfMC(t,:) = xparticlePred_pfukfMC(t,outIndex); % Keep particles with
                                              % resampled indices.
  Pparticle_pfukfMC(t,:) = PparticlePred_pfukfMC(t,outIndex); 
  
  previousXResukfMC(t,:) = previousXukfMC(t,outIndex);  % Resample particles
                                                        % at t-1.
  previousPResukfMC(t,:) = previousPukfMC(t,outIndex);  % Resample particles
                                                        % at t-1.
   
  % METROPOLIS-HASTINGS STEP:
  % ========================
  u=rand(N,1); 
  accepted=0;
  rejected=0;
  for i=1:N,   
     % Call Unscented Kalman Filter
    [muProp,PProp]=ukf(previousXResukfMC(t,i),previousPResukfMC(t,i),[],Q_pfukf,'ukf_ffun',y(t),R_pfukf,'ukf_hfun',t,alpha,beta,kappa);    
    xparticleProp = muProp + sqrtm(PProp)*randn(1,1);
    PparticleProp = PProp;   
    mProp = feval('hfun',xparticleProp,t);        
    likProp = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-mProp)^(2)))+1e-99;
    priorProp = ((xparticleProp-previousXResukfMC(t,i))^(g1-1)) ...
		 * exp(-g2*(xparticleProp-previousXResukfMC(t,i)));
    proposalProp = inv(sqrt(PparticleProp)) * ...
	       exp(-0.5*inv(PparticleProp) *( ...
					      (xparticleProp-muProp)^(2)));
    m = feval('hfun',xparticle_pfukfMC(t,i),t);        
    lik = inv(sqrt(sigma)) * exp(-0.5*inv(sigma)*((y(t)-m)^(2)))+1e-99;
    prior = ((xparticle_pfukfMC(t,i)-previousXResukfMC(t,i))^(g1-1)) ...
		 * exp(-g2*(xparticle_pfukfMC(t,i)-previousXResukfMC(t,i)));
    proposal = inv(sqrt(Pparticle_pfukfMC(t,i))) * ...
	       exp(-0.5*inv(Pparticle_pfukfMC(t,i)) *((xparticle_pfukfMC(t,i)-muProp)^(2)));
    ratio = (likProp*priorProp*proposal)/(lik*prior*proposalProp);
    acceptance = min(1,ratio);
    if u(i,1) <= acceptance 
      xparticle_pfukfMC(t,i) = xparticleProp;
      Pparticle_pfukfMC(t,i) = PparticleProp;
      accepted=accepted+1;
    else
      xparticle_pfukfMC(t,i) = xparticle_pfukfMC(t,i); 
      Pparticle_pfukfMC(t,i) = Pparticle_pfukfMC(t,i);  
      rejected=rejected+1;
    end;
  end;   % End of MCMC loop. 
end;   % End of t loop.

time_pfukfMC(j) = toc;



%%%%%%%%%%%%%%%%%%%%%  PLOT THE RESULTS  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  ================  %%%%%%%%%%%%%%%%%%%%%

figure(1)
clf;
p0=plot(1:T,y,'k+','lineWidth',2); hold on;
%p2=plot(1:T,mu_ekf,'r:','lineWidth',2); hold on;
%p3=plot(1:T,mu_ukf,'b:','lineWidth',2);
p4=plot(1:T,mean(xparticle_pf(:,:)'),'g','lineWidth',2);
p5=plot(1:T,mean(xparticle_pfekf(:,:)'),'r','lineWidth',2);
p6=plot(1:T,mean(xparticle_pfukf(:,:)'),'b','lineWidth',2); 
p1=plot(1:T,x,'k:o','lineWidth',2); hold off;
%legend([p1 p2 p3 p4 p5 p6],'True x','EKF estimate','UKF estimate','PF estimate','PF-EKF estimate','PF-UKF estimate');
legend([p0 p1 p4 p5 p6],'Noisy observations','True x','PF estimate','PF-EKF estimate','PF-UKF estimate');
xlabel('Time','fontsize',15)
zoom on;
title('Filter estimates (posterior means) vs. True state','fontsize',15)

figure(2)
clf
subplot(211);
semilogy(1:T,P_ekf,'r--',1:T,P_ukf,'b','lineWidth',2);
legend('EKF','UKF');
title('Estimates of state covariance','fontsize',14);
xlabel('time','fontsize',12);
ylabel('var(x)','fontsize',12);
zoom on;

if (1),
figure(3)
clf;
% Plot predictive distribution of y:
subplot(231);
domain = zeros(T,1);
range = zeros(T,1);
thex=[-3:.1:15];
hold on
ylabel('Time (t)','fontsize',15)
xlabel('y_t','fontsize',15)
zlabel('p(y_t|y_{t-1})','fontsize',15)
title('Particle Filter','fontsize',15);
%v=[0 1];
%caxis(v);
for t=6:5:T,
  [range,domain]=hist(yPred_pf(t,:),thex);
  waterfall(domain,t,range/sum(range));
end;
view(-30,80);
rotate3d on;
a=get(gca);
set(gca,'ygrid','off');
% Plot posterior distribution of x:
subplot(234);
domain = zeros(T,1);
range = zeros(T,1);
thex=[0:.1:10];
hold on
ylabel('Time (t)','fontsize',15)
xlabel('x_t','fontsize',15)
zlabel('p(x_t|y_t)','fontsize',15)
%v=[0 1];
%caxis(v);
for t=6:5:T,
  [range,domain]=hist(xparticle_pf(t,:),thex);
  waterfall(domain,t,range/sum(range));
end;
view(-30,80);
rotate3d on;
a=get(gca);
set(gca,'ygrid','off');

% Plot predictive distribution of y:
subplot(232);
domain = zeros(T,1);
range = zeros(T,1);
thex=[-3:.1:15];
hold on
ylabel('Time (t)','fontsize',15)
xlabel('y_t','fontsize',15)
zlabel('p(y_t|y_{t-1})','fontsize',15)
title('Particle Filter (EKF proposal)','fontsize',15);
%v=[0 1];
%caxis(v);
for t=6:5:T,
  [range,domain]=hist(yPred_pfekf(t,:),thex);
  waterfall(domain,t,range/sum(range));
end;
view(-30,80);
rotate3d on;
a=get(gca);
set(gca,'ygrid','off');
% Plot posterior distribution of x:
subplot(235);
domain = zeros(T,1);
range = zeros(T,1);
thex=[0:.1:10];
hold on
ylabel('Time (t)','fontsize',15)
xlabel('x_t','fontsize',15)
zlabel('p(x_t|y_t)','fontsize',15)
%v=[0 1];
%caxis(v);
for t=6:5:T,
  [range,domain]=hist(xparticle_pfekf(t,:),thex);
  waterfall(domain,t,range/sum(range));
end;
view(-30,80);
rotate3d on;
a=get(gca);
set(gca,'ygrid','off');

% Plot predictive distribution of y:
subplot(233);
domain = zeros(T,1);
range = zeros(T,1);
thex=[-3:.1:15];
hold on
ylabel('Time (t)','fontsize',15)
xlabel('y_t','fontsize',15)
zlabel('p(y_t|y_{t-1})','fontsize',15)
title('Particle Filter (UKF proposal)','fontsize',15);
%v=[0 1];
%caxis(v);
for t=6:5:T,
  [range,domain]=hist(yPred_pfukf(t,:),thex);
  waterfall(domain,t,range/sum(range));
end;
view(-30,80);
rotate3d on;
a=get(gca);
set(gca,'ygrid','off');
% Plot posterior distribution of x:
subplot(236);
domain = zeros(T,1);
range = zeros(T,1);
thex=[0:.1:10];
hold on
ylabel('Time (t)','fontsize',15)
xlabel('x_t','fontsize',15)
zlabel('p(x_t|y_t)','fontsize',15)
%v=[0 1];
%caxis(v);
for t=6:5:T,
  [range,domain]=hist(xparticle_pfukf(t,:),thex);
  waterfall(domain,t,range/sum(range));
end;
view(-30,80);
rotate3d on;
a=get(gca);
set(gca,'ygrid','off');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- CALCULATE PERFORMANCE --%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmsError_ekf(j)     = sqrt(inv(T)*sum((x-mu_ekf).^(2)));
rmsError_ukf(j)     = sqrt(inv(T)*sum((x-mu_ukf).^(2)));
rmsError_pf(j)      = sqrt(inv(T)*sum((x'-mean(xparticle_pf')).^(2)));
rmsError_pfMC(j)    = sqrt(inv(T)*sum((x'-mean(xparticle_pfMC')).^(2)));
rmsError_pfekf(j)   = sqrt(inv(T)*sum((x'-mean(xparticle_pfekf')).^(2)));
rmsError_pfekfMC(j) = sqrt(inv(T)*sum((x'-mean(xparticle_pfekfMC')).^(2)));
rmsError_pfukf(j)   = sqrt(inv(T)*sum((x'-mean(xparticle_pfukf')).^(2)));
rmsError_pfukfMC(j) = sqrt(inv(T)*sum((x'-mean(xparticle_pfukfMC')).^(2)));

disp(' ');
disp('Root mean square (RMS) errors');
disp('-----------------------------');
disp(' ');
disp(['EKF          = ' num2str(rmsError_ekf(j))]);
disp(['UKF          = ' num2str(rmsError_ukf(j))]);
disp(['PF           = ' num2str(rmsError_pf(j))]);
disp(['PF-MCMC      = ' num2str(rmsError_pfMC(j))]);
disp(['PF-EKF       = ' num2str(rmsError_pfekf(j))]);
disp(['PF-EKF-MCMC  = ' num2str(rmsError_pfekfMC(j))]);
disp(['PF-UKF       = ' num2str(rmsError_pfukf(j))]);
disp(['PF-UKF-MCMC  = ' num2str(rmsError_pfukfMC(j))]);

disp(' ');
disp(' ');
disp('Execution time  (seconds)');
disp('-------------------------');
disp(' ');
disp(['PF           = ' num2str(time_pf(j))]);
disp(['PF-MCMC      = ' num2str(time_pfMC(j))]);
disp(['PF-EKF       = ' num2str(time_pfekf(j))]);
disp(['PF-EKF-MCMC  = ' num2str(time_pfekfMC(j))]);
disp(['PF-UKF       = ' num2str(time_pfukf(j))]);
disp(['PF-UKF-MCMC  = ' num2str(time_pfukfMC(j))]);
disp(' ');

drawnow;

%*************************************************************************

end    % Main loop (for j...)

% calculate mean of RMSE errors
mean_RMSE_ekf     = mean(rmsError_ekf);
mean_RMSE_ukf     = mean(rmsError_ukf);
mean_RMSE_pf      = mean(rmsError_pf);
mean_RMSE_pfMC    = mean(rmsError_pfMC);
mean_RMSE_pfekf   = mean(rmsError_pfekf);
mean_RMSE_pfekfMC = mean(rmsError_pfekfMC);
mean_RMSE_pfukf   = mean(rmsError_pfukf);
mean_RMSE_pfukfMC = mean(rmsError_pfukfMC);

% calculate variance of RMSE errors
var_RMSE_ekf     = var(rmsError_ekf);
var_RMSE_ukf     = var(rmsError_ukf);
var_RMSE_pf      = var(rmsError_pf);
var_RMSE_pfMC    = var(rmsError_pfMC);
var_RMSE_pfekf   = var(rmsError_pfekf);
var_RMSE_pfekfMC = var(rmsError_pfekfMC);
var_RMSE_pfukf   = var(rmsError_pfukf);
var_RMSE_pfukfMC = var(rmsError_pfukfMC);

% calculate mean of execution time
mean_time_pf      = mean(time_pf);
mean_time_pfMC    = mean(time_pfMC);
mean_time_pfekf   = mean(time_pfekf);
mean_time_pfekfMC = mean(time_pfekfMC);
mean_time_pfukf   = mean(time_pfukf);
mean_time_pfukfMC = mean(time_pfukfMC);

% display final results

disp(' ');
disp(' ');
disp('************* FINAL RESULTS *****************');
disp(' ');
disp('RMSE : mean and variance');
disp('---------');
disp(' ');
disp(['EKF          = ' num2str(mean_RMSE_ekf) ' (' num2str(var_RMSE_ekf) ')']);
disp(['UKF          = ' num2str(mean_RMSE_ukf) ' (' num2str(var_RMSE_ukf) ')']);
disp(['PF           = ' num2str(mean_RMSE_pf) ' (' num2str(var_RMSE_pf) ')']);
disp(['PF-MCMC      = ' num2str(mean_RMSE_pfMC) ' (' num2str(var_RMSE_pfMC) ')']);
disp(['PF-EKF       = ' num2str(mean_RMSE_pfekf) ' (' num2str(var_RMSE_pfekf) ')']);
disp(['PF-EKF-MCMC  = ' num2str(mean_RMSE_pfekfMC) ' (' num2str(var_RMSE_pfekfMC) ')']);
disp(['PF-UKF       = ' num2str(mean_RMSE_pfukf) ' (' num2str(var_RMSE_pfukf) ')']);
disp(['PF-UKF-MCMC  = ' num2str(mean_RMSE_pfukfMC) ' (' num2str(var_RMSE_pfukfMC) ')']);

disp(' ');
disp(' ');
disp('Execution time  (seconds)');
disp('-------------------------');
disp(' ');
disp(['PF           = ' num2str(mean_time_pf)]);
disp(['PF-MCMC      = ' num2str(mean_time_pfMC)]);
disp(['PF-EKF       = ' num2str(mean_time_pfekf)]);
disp(['PF-EKF-MCMC  = ' num2str(mean_time_pfekfMC)]);
disp(['PF-UKF       = ' num2str(mean_time_pfukf)]);
disp(['PF-UKF-MCMC  = ' num2str(mean_time_pfukfMC)]);
disp(' ');

%*************************************************************************

break;

% This is an alternative way of plotting the 3D stuff:
% Somewhere in between lies the best way!
figure(3)
clf;
support=[-1:.1:2];
NN=50;
extPlot=zeros(10*61,1);
for t=6:5:T,
  [r,d]=hist(yPred_pf(t,:),support);
  r=r/sum(r);
  for i=1:1:61,
    for j=1:1:NN,
      extPlot(NN*i-NN+1:i*NN) = r(i);
    end;
  end;
  d= linspace(-5,25,length(extPlot));
  plot3(d,t*ones(size(extPlot)),extPlot,'r')
  hold on;
end;



















