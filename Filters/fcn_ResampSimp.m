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
% function outIndex = fcn_ResampSimp(w, N)
% Draw a total of N samples with probabilities proportional to the
% weight vector w, using sorting.
% w            : normalized weight vector (sum to one)
% N (optional) : total number of samples; default to length(w)
% outIndex     : each element is an index into w, or, the "parent" of
%                the sample. Therefore if {X, w} is the original
%                particles-weights pair, then {X(outIndex), 1/N}
%                will be the resampled pair.
%
% See also     : fcn_ResampSys, fcn_ResampResid
function outIndex = fcn_ResampSimp(w, N)

eps = 1e-12; % small but not too small

len = length(w);
F = cumsum(w);
if abs(F(end) - 1) > eps
    error('the weight vector should be normalized.');
end;

switch nargin
    case 1,
        N = len;
    case 2,
    otherwise,
        error('wrong number of arguments');
end;

s = sort(rand(1, N)); % faster than cumprod(rand(1,N).^(1./[N:-1:1]))

outIndex = zeros(1, N);
j = 1;
for i = 1:N
    while F(j) < s(i)
        j = j + 1;
    end;
    outIndex(i) = j;
end;

end
