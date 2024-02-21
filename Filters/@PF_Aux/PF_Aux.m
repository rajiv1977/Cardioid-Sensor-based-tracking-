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
% function fltr = PF_Aux(sys, init_distr, opts)
% PF_Aux: class constructor for an Auxiliary  Particle Filter, with a simple
% proposal (the state transition) but modified resampling.
%
% sys       : a SigObsSys object
% init_distr: initial distribution of the particles
% opts      : a structure array of options, obtained from calling setPFOptions
% fltr      : the contructed filter
function fltr = PF_Aux(sys, init_distr, opts)

if nargin == 0
    fltr.f = []; fltr.h = []; fltr.w_d = []; fltr.v_d = [];
    fltr.N = []; fltr.phi = [];
    fltr.p = []; fltr.w = [];
    fltr = class(fltr,'PF_Aux');
elseif isa(sys, 'PF_Aux')
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
    [x, y, f, h, w_d, v_d, phi] = get(sys);
    fltr.f = f; % state equation
    fltr.h = h; % measurement equation
    fltr.w_d = w_d; % distribution of state noise
    fltr.v_d = v_d; % distribution of measurement noise
    fltr.phi = phi;
    try
        fltr.N = opts.NumPart;
    catch
        error(['wrong option structuure:' lasterr]);
    end;
    
    try
        fltr.p = drawSamples(init_distr, fltr.N); % particles
        fltr.w = ones(1, fltr.N) / fltr.N;        % weights
    catch
        error('wrong initial particle distribution');
    end;
    fltr = class(fltr, 'PF_Aux');
end

end