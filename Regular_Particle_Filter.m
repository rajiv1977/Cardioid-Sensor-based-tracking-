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
clc;
clear;
clear classes
tic;
current_dir_name = cd;
addpath([current_dir_name '/Filters']);

% The dimension variables were used by PFL for error checking
dimX = 3;
dimY = 3;
T               = 0.0001;     % Sampling Time
Nsim            = 1024;     % Number of scans
num_particles   = 600;      % Number of particles

% initial condition of the system; NOT used to initialize the filter
initX = [2 30*(pi/180) 120];

% function handle for state transition
f = @(x) [x(1); x(2); x(3)];

Phase_ang = (pi/180)*randi([-45 45],1,1);

% function handle for observation
h = @(x, ti, phi) [  2*x(1) * cos( (x(3)*ti) + phi);
    2*x(1) * cos( (x(3)*ti) + phi) * cos(x(2));
    2*x(1) * cos( (x(3)*ti) + phi) * sin(x(2))];

meanW = [0 0 0];
p_noise_amp = 0.01;
p_noise_ber = 2.5*(pi/180);
p_noise_frq = 1;
p_noise_var     = [p_noise_amp p_noise_ber p_noise_frq];
covW = eye(3);
covW(1,1) = p_noise_var(1);
covW(2,2) = p_noise_var(2);
covW(3,3) = p_noise_var(3);

q=0.05;
covW = covW*q;

meta = 1/2;
meas_var = 0.05;

R  = [ 1    0       0;
    0    meta    0;
    0    0       meta] * meas_var;
P_ini           = R;

%covW = eye(3);
% constructing a Gaussian distribution
wDistr =  GaussianDistr(meanW, covW);

meanV = [0 0 0];
covV = eye(3);
covV(1,1) = sqrt(R(1,1));
covV(2,2) = sqrt(R(2,2));
covV(3,3) = sqrt(R(3,3));

vDistr =  GaussianDistr(meanV, covV);

% Constructing the system postulated by the filter.
theSys = SigObsSys(f, h, wDistr, vDistr, initX, T, Phase_ang);

initMean = initX + (sqrtm(P_ini)*randn(3,1))';
initCov = P_ini;

% constructing the initial distribution postulated by the filter
initDistr = GaussianDistr(initMean, initCov);

% constructing the options
opts = setPFOptions('NumPart', num_particles, 'ResampAlgo', 'fcn_ResampResid', 'ResampPeriod',2);

% constructing the particle filter
theFltr = PF_Simple(theSys, initDistr, opts);
 
% Here is an example of a simulation.
% Assuming that the real system from which data is collected is given as:
realSys = theSys;

% % We can reset the initial condition
% realSys = set(realSys, rand(dimX, 1));

% Now we simulate the signal and observation trajectories.
X = zeros(dimX, Nsim); Y = zeros(dimY, Nsim); Xhat = X;
ti = 0;
b_hat(1) = 0;
for i = 2:Nsim
    ti = ti + T;
    Time_index(i) = ti;
    realSys = update(realSys, ti);
    [X(:, i), Y(:, i)] = get(realSys);
    ang = arctangent_estimator(Y,T);
    b_hat = [b_hat; ang];
end;

% Now we simulate the filtering.
ti = 0;
ErR_arcT = zeros(Nsim,1);
for i=2:Nsim
    i
    ti = ti + T;
    theFltr = update(theFltr, Y(:, i), ti);
    [Xhat(:, i), scan(i).Phat, scan(i).pt, scan(i).wt] = get(theFltr);
    ErR =  X(:,i) - Xhat(:, i);
    nees(i) = ErR' * inv(scan(i).Phat) * ErR;
    rmse_ekf_bea(i) = sqrt((ErR(2))^2);
    rmse_act_bea(i) = sqrt((X(2,i) - b_hat(i))^2);
end;
 
figure(4)
Tm = Time_index;
Tm(1) = [];
Xtruth = X(2,:);
Xtruth(1) = [];
b_hat(1) = [];
plot(Tm,Xtruth); hold on; plot(Time_index,Xhat(2,:), 'r-.'); hold on; plot(Tm,b_hat, 'g-.');
legend('True bearing', 'EKF - Estimated bearing', 'arctangent estimator - Estimated bearing');
xlabel('Time_indexe (s)');
ylabel('Bearing (rad) ');

figure(5)
plot(Time_index,X(1,:)); hold on; plot(Time_index,Xhat(1,:), 'r-.');
legend('True amplitude','Estimated amplitude');
xlabel('Time_indexe (s)');
ylabel('Amplitude ');

figure(6)
plot(Time_index,X(2,:)); hold on; plot(Time_index,Xhat(2,:), 'r-.');
legend('True bearing','Estimated bearing');
xlabel('Time (s)');
ylabel('Bearing (rad)');

figure(7)
plot(Time_index,X(3,:)); hold on; plot(Time_index,Xhat(3,:), 'r-.');
legend('True frequency','Estimated frequency');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

figure(8);
plot(Time_index,nees,'--r','LineWidth',1,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',1)
hold on;
lower_bound = ones(1,length(nees))*0.352;
upper_bound = ones(1,length(nees))*7.81;
plot(Time_index,lower_bound,'--k','LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',1)
hold on;
plot(Time_index,upper_bound,'--k','LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',1)
legend('NEES', '95% Confidence region');
xlabel('Time (s)');
ylabel('NEES');

figure(9);
plot(Time_index,rmse_ekf_bea,'--r','LineWidth',1,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',1)
hold on;
plot(Time_index,rmse_act_bea,'--b','LineWidth',1,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',1)
hold off;
legend('EKF - RMSE', 'Arctangent RMSE');
xlabel('Time (s)');
ylabel('RMSE (bearing)');


toc;
% X
% Xhat