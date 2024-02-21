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
% function fltr = PF_EKF(sys, init_distr, opts)
% PF_EKF: class constructor for a Particle Filter with a proposal
% determined by an EKF-type procedure
%
% sys       : a SigObsSys object, with Gaussian noises, and
%             Jacobian for the observation
% init_distr: initial distribution of the particles
% opts      : a structure array of options, obtained from calling setPFOptions
% fltr      : the contructed filter
function fltr = PF_EKF(sys, init_distr, opts)

if nargin == 0
    fltr.f = []; fltr.h = []; fltr.w_d = []; fltr.v_d = [];
    fltr.invQ = []; fltr.invR = [];
    fltr.df = []; fltr.dh = [];
    fltr.N = []; fltr.T = []; fltr.resamp = []; fltr.thresh = [];
    fltr.p = []; fltr.w = []; fltr.counter = [];
    fltr.phi = [];
    fltr = class(fltr,'PF_EKF');
elseif isa(sys, 'PF_EKF')
    fltr = sys;
elseif nargin ~= 3
    error('wrong number of arguments');
else
    if ~isa(sys, 'SigObsSys')
        error('wrong first argument: must be a SigObsSys object');
    end;
    if ~setPFOptions(opts)
        error(['wrong second argument: must be a valid option' ...
            ' structure']);
    end;
    [x, y, f, h, w_d, v_d, phi, df, dh] = get(sys);
    
    fltr.f = f; % state equation
    fltr.h = h; % measurement equation
    fltr.w_d = w_d; % distribution of state noise
    fltr.v_d = v_d; % distribution of measurement noise
    fltr.invQ = inv(cov(w_d)); % constant for EKF proposal
    % calculation
    fltr.invR = inv(cov(v_d)); % ditto
    fltr.df = df; % Jacobian of state function, can be []
    fltr.dh = dh; % Jacobian of observation function, needed for EKF
    fltr.phi = phi;
    % proposal
    
    try
        fltr.N = opts.NumPart;
        fltr.T = opts.ResampPeriod;
        
        if strmatch(opts.ResampAlgo, 'none', 'exact')
            fltr.r_fun = [];
        else
            eval(['fltr.r_fun = @' opts.ResampAlgo ';']);
        end;
        
        if strmatch(opts.ResampAlgo, 'fcn_ResampBran', 'exact')
            fltr.thresh = opts.BranchThresh;
        else fltr.thresh = 1;
        end;
    catch
        error(['wrong option structuure:' lasterr]);
    end;
    
    try
        fltr.p = drawSamples(init_distr, fltr.N); % particles
        fltr.w = ones(1, fltr.N) / fltr.N;        % weights
    catch
        error('wrong initial particle distribution');
    end;
    fltr.counter = 0; % internal counter for determining when to resample
    fltr = class(fltr, 'PF_EKF');
end

end
