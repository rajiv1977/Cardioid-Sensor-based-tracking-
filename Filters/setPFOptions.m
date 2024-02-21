%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Cardioid Sensor based tracking                                  %
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
% function opts = setPFOptions(varargin)
% Create an options structure for initializing Particle Filter
% classes
%   varargin: empty, or parameter-value pairs, i.e, 'param1', 'value1',
%             'param2', 'value2', etc. Parameter name is case-insensitive.
%   opts:     the options structure
% Called without input and output arguments: print a list of
% options
% Called without input argument but with output argument: return
% the default
% Called with parameter-value pairs: replace the corresponding
% default values
% Note: Each Particle Filter class may use only a subset of the
%   parameters created by this function. 
function opts = setPFOptions(varargin)

if nargin == 0 && nargout == 0 % print help
  fprintf(' NumPart      : [ positive integer : {1000} ]\n');
  fprintf(' ResampPeriod : [ positive integer : {10} ]\n');
  fprintf([' ResampAlgo   : [ ''none'' | ''fcn_ResampSimp'' | ''fcn_ResampResid''' ...
	   ' | ''fcn_ResampBran''| {''fcn_ResampSys''} ]\n']);
  fprintf(' BranchThresh : [ positive scalar between 0 and 1 : {0.5} ]\n');
  fprintf(' RegularAlgo  : [{''post''}, | ''pre'' | ''mix'']\n');
  fprintf(' RegularWidth : [ positive scalar : {0.1}]\n');
  fprintf(' RegularWhitening : [ boolean : {0} | 1]\n');
  return;
end;


params = {'NumPart', 'ResampPeriod', 'ResampAlgo', 'BranchThresh', ...
	  'RegularAlgo', 'RegularWidth', 'RegularWhitening'};
values = {1000, 10, 'fcn_ResampSys', 0.5, 'post', 0.1, 0};
tmp = [params; values];
opts = struct(tmp{:});

if nargin == 0 && nargout == 1 % return default
  return;
end;

if nargin == 1 % check validity of the entire structure. This is a
               % hack so that the "is_valid(p,v)" subfunction can
               % be reused, instead of being pulled out. 
  try 
    names = fieldnames(varargin{1});
    opts = 1;
    for i = 1:length(names)
      if ~is_valid(lower(names{i}), getfield(varargin{1}, names{i}))
	opts = 0;
	break;
      end;
    end;
  catch
    error('wrong option structure');
  end;
  return;
end;  

if rem(nargin, 2) ~=0 % set options
  error('arguments must be in parameter-value pairs');
end;

lowparams = lower(params);
i = 1;
while i <= nargin
  p = varargin{i};
  v = varargin{i+1};
  if ~ischar(p)
    error('argument %d: must be string', i);
  end;
  lp = lower(strtrim(p));
  
  j = strmatch(lp, lowparams, 'exact');
  if isempty(j)
    error('Argument %d: wrong name.', i);
  elseif is_valid(lowparams{j}, v)
    opts.(params{j}) = v;
  else 
    error('Argument %d: wrong value.', i+1);
  end;
  i = i + 2;
end;

function result = is_valid(p, v)

result = 0;
switch p
 case {'numpart', 'resampperiod'},
  if isreal(v) && isscalar(v) && v > 0 && abs(v - round(v)) < eps
    result = 1;
  end;
 case {'resampalgo'},
  if ischar(v) && ~isempty(strmatch(v, {'none', 'fcn_ResampSimp', ...
		    'fcn_ResampResid', 'fcn_ResampBran', 'fcn_ResampSys'}, 'exact'))
    result = 1;
  end;
 case {'branchthresh'},
  if isreal(v) && isscalar(v) && v > 0 && v < 1
    result = 1;
  end;
 case {'regularalgo'},
  if ischar(v) && ~isempty(strmatch(v, {'pre', 'post', 'mix'}, ...
				   'exact'))
    result = 1;
  end;
 case {'regularwidth'},
  if isreal(v) && isscalar(v) && v > 0
    result = 1;
  end;
 case {'regularwhitening'},
  if isreal(v) && isscalar(v)
    result = 1;
  end;
end;

