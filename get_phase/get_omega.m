function outFreq=get_omega(PHIc,sr,m,g)
% get median of input signal derivative 
% input:
%       signal: input signal
%       sr: sampling rate

%       m: length of the Savitzky-Golay differentiator
%       g: coefficients of the differentiatior
% output:
%       allFreqs: frequency values
%       Ts: vector of time steps

% Leonardo Lancia 14/10/2023
% mailto: leonardo.lancia@cnrs.fr

if nargin<3 || isempty(m)% initialize m
    m=5;
end
if nargin<4 || isempty(g)% initialize g
    g = SG_calc(m,3, 1)'; %compute coeffs for the sgolay differentiator
end


p=1; % index of the relevant S.G. coefficients
h=1/sr; % sampling period
M = (m+1)/2;      % output point , M = (m+1)/2; % middle point for odd m

myKernel=factorial(p)/(-h)^p * g(:,p+1); % cernel of the differentiator


PHIu=unwrap(PHIc); % unwrap phase

dPHI = filter(myKernel,1,PHIu); % differentiation
%dPHI=-sample_advance(dPHI,m/2); % correct for filter delay
dPHI=-[dPHI(m/2+1:end)]./sr;% correct for filter delay without padding with zeros
outFreq=median(dPHI./(2*pi));%get median frequency


