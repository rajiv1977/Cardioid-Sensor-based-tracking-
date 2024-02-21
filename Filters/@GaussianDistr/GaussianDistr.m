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
% function gDistr = GaussianDistr(theMean, theCov)
% GaussianDistr: Gaussian distribution class constructor.
% gDistr = GaussianDistr(theMean, theCov)
% creates an object that represents a Gaussian (i.e., Normal)
% distribution with theMean as the mean vector and theCov as the
% covariance matrix. The covariance matrix must be positive definite in
% this implementation.
function gDistr = GaussianDistr(theMean, theCov)
if nargin == 0
    gDistr.mean = [];
    gDistr.n = [];
    gDistr.cov = [];
    gDistr.const = [];
    gDistr.invCov = [];
    gDistr.UsqrtT = [];
    gDistr = class(gDistr,'GaussianDistr');
elseif isa(theMean,'GaussianDistr')
    gDistr = theMean;
elseif nargin == 2
    gDistr.mean = theMean(:);
    gDistr.n = length(gDistr.mean);
    [nrow, ncol] = size(theCov);
    if nrow ~= ncol || nrow ~= gDistr.n
        error('wrong size for mean or covariance');
    end;
    tmp = theCov';
    if any(theCov(:) ~= tmp(:))
        error('covariance matrix should be symmetric');
    end;
    [U, T] = schur(theCov);
    if any(diag(T) <= eps)
        error('covariance matrix should be positive definite');
    end;
    gDistr.cov = theCov;
    % constant used in the density calculation
    gDistr.const = 1 / sqrt((2 * pi)^gDistr.n * det(theCov));
    % pre-computed inverse for density calculation
    gDistr.invCov = inv(theCov);
    % pre-computed matrices for drawing samples
    gDistr.UsqrtT = U * (T .^ 0.5);
    gDistr = class(gDistr, 'GaussianDistr');
else
    error('wrong number of arguments');
end

end
