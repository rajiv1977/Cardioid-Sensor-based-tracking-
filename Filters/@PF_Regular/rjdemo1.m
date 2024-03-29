%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Cardioid Sensor based tracking                                  %
%                            Regularized Particle Filter                                    %
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
% PURPOSE : To approximate a noisy nolinear function with RBFs, where the number
%           of parameters and parameter values are estimated via reversible jump
%           Markov Chain Monte Carlo (MCMC) simulation.
clear;
echo off;

% INITIALISATION AND PARAMETERS:
% =============================

N = 100;                       % Number of time steps.
t = 1:1:N;                    % Time.
chainLength = 2000;           % Length of the Markov chain simulation.
burnIn = 1000;                % Burn In period.
bFunction = 'rjGaussian'      % Type of basis function.
par.doPlot = 0;               % 1 plot. 0 don't plot.
par.a = 2;                    % Hyperparameter for delta.  
par.b = 10;                   % Hyperparameter for delta.
par.e1 = 0.0001;              % Hyperparameter for nabla.    
par.e2 = 0.0001;              % Hyperparameter for nabla.   
par.v = 0;                    % Hyperparameter for sigma    
par.gamma = 0;                % Hyperparameter for sigma. 
par.kMax = 50;                % Maximum number of basis.
par.arbC = 0.25;              % Constant for birth and death moves.
par.merge = .1;               % Split-Merge parameter.
par.Lambda = .5;              % Hybrid Metropolis decision parameter.
par.sRW = .001;               % Variance of noise in the random walk.
par.walkPer = 0.1;            % Percentange of random walk interval. 

% GENERATE THE DATA:
% =================
noiseVar = 0.5;
x = 4*rand(N,1)-2;                    % Input data - uniform in [-2,2].
u = randn(N,1);
noise = sqrt(noiseVar)*u;             % Measurement noise
varianceN=var(noise)
y = x + 2*exp(-16*(x.^(2))) + 2*exp(-16*((x-.7).^(2)))   + noise;  % Output data.
x=(x+2)/4;                            % Rescaling to [0,1].
ynn = y-noise;
xv = 4*rand(N,1)-2;                    % Input data - uniform in [-2,2].
uv = randn(N,1);
noisev = sqrt(noiseVar)*uv;   
yv = xv + 2*exp(-16*(xv.^(2))) + 2*exp(-16*((xv-.7).^(2)))   + noisev;  % Output data.
xv=(xv+2)/4;               
yvnn = yv-noisev;

figure(1)
subplot(211)
plot(x,y,'b+');
ylabel('Output data','fontsize',15);
xlabel('Input data','fontsize',15);
%axis([0 1 -3 3]);
subplot(212)
plot(noise)
ylabel('Measurement noise','fontsize',15);
xlabel('Time','fontsize',15);

fprintf('\n')
fprintf('Press a key to continue')
fprintf('\n')
pause;

% PERFORM REVERSE JUMP MCMC WITH RADIAL BASIS:
% ===========================================
[k,mu,alpha,sigma,nabla,delta,yp,ypv,post] = rjnn(x,y,chainLength,N,bFunction,par,xv,yv);

% COMPUTE CENTROID, MAP AND VARIANCE ESTIMATES:
% ============================================

[l,m]=size(mu{1});
[Nv,d]=size(xv);
l=chainLength-burnIn;
muvec=zeros(l,m);
alphavec=zeros(m+d+1,l);
ypred = zeros(N,l+1);
ypredv = zeros(Nv,l+1);
for i=1:N;
  ypred(i,:) = yp(i,1,burnIn:chainLength);
end;
for i=1:Nv;
  ypredv(i,:) = ypv(i,1,burnIn:chainLength);
end;
ypred = mean(ypred');
ypredv = mean(ypredv');
fevTrain =(y-ypred')'*(y-ypred')*inv((y-mean(y)*ones(size(y)))'*(y-mean(y)*ones(size(y))))
fevTest = (yv-ypredv')'*(yv-ypredv')*inv((yv-mean(yv)*ones(size(yv)))'*(yv-mean(yv)*ones(size(yv))))

% PLOTS:
% =====
figure(1)
clf;
[xv,i]=sort(xv);
yvnn=yvnn(i);
ypredv=ypredv(i);
yv=yv(i);
[x,i]=sort(x);
ynn=ynn(i);
ypred=ypred(i);
y=y(i);
clf;
subplot(211)
plot(x,ynn,'m--',x,y,'b+',x,ypred,'r')
ylabel('Train output','fontsize',15)
xlabel('Train input','fontsize',15)
subplot(212)
plot(xv,yvnn,'m--',xv,yv,'b+',xv,ypredv,'r')
ylabel('Test output','fontsize',15)
xlabel('Test input','fontsize',15)
legend('True function','Test data','Prediction');

% COMPUTE THE MOST LIKELY MODES:
% =============================
pInt=2;
support=[1:1:4];
probk=zeros((chainLength)/pInt,length(support));
for p=pInt:pInt:chainLength,  
  [probk(p/pInt,:),kmodes]=hist(k(1:p),support);
  probk(p/pInt,:)=probk(p/pInt,:)/p;
end;
figure(2)
clf;
plot(pInt:pInt:chainLength,probk(:,1),'k--',pInt:pInt:chainLength,probk(:,2),'k:',pInt:pInt:chainLength,probk(:,3),'k',pInt:pInt:chainLength,probk(:,4),'k-.');
xlabel('Chain length','fontsize',15)
ylabel('p(k|y)','fontsize',15)
legend('1','2','3','4')
modes = probk(chainLength/2,:);

% HISTOGRAMS:
% ==========
figure(3)
clf;
subplot(321)
hist(delta(burnIn:chainLength),80)
ylabel('Regularisation parameter','fontsize',15);
subplot(322)
plot(delta)
ylabel('Regularisation parameter','fontsize',15);
subplot(323)
hist(sigma(burnIn:chainLength),80)
ylabel('Noise variance','fontsize',15);
subplot(324)
plot(sigma)
ylabel('Noise variance','fontsize',15);
subplot(325)
hist(nabla(burnIn:chainLength),80)
ylabel('Poisson parameter','fontsize',15);
subplot(326)
plot(nabla)
ylabel('Poisson parameter','fontsize',15);















