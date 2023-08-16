function [PHI,PHIc]=getPHI(signal,sr,m,n,centering,handleNegPHI)

%input: 
% signal: input signal
% sr: sampling rate
% m: number of filtered points for Savitzky-Golay differentiator 
% n: polinomial order of the differentiator
%centering: parameter handling centering strategy. If EMD, centering corresponds to extracting the first EMD component. 
% handleNegPHI: parameter handling negative frequency correction. Can be interp or a scalar, in
%   the former case the values are interpolated otherwise they are replaced
%   with the scalar
%output:
% PHI: instantaneous phase
% PHIc: instantaneous phase corrected for neg frequencies

% example1 (don't use EMD, interpolate phase vals corresponding to neg freqs
% phase difference (freq) is computed with a polinomial order = 3 and on 5 points at a time:
%[PHI,PHIc]=getPHI(sigInf,sr,5,3,'','interp');
%example 2 (equivalent to example 1 because default behaviour):
%[PHI,PHIc]=getPHI(sigInf,sr,5,3)
% example3 (use EMD, everything else being equal to other examples):
%[PHI,PHIc]=getPHI(sigInf,sr,5,3,'EMD','interp');


if nargin<3
    m=5;% number of filtered points (FIR order)
end
if nargin<4
    n=3;% approximation polynomial order
end
if nargin<5
    centering=''
end
if nargin<6
    handleNegPHI='interp';
end

M = (m+1)/2;      % output point , M = (m+1)/2; % middle point for odd m
h = 1/sr; %sample period
SG_coef = SG_calc(m,n,h); %compute coeffs for the sgolay differentiator
origLen=length(signal); %original length of the signal
tmpLen=2^nextpow2(origLen);%set optimal signal length to go through the Hilbert transform

if strcmp(centering,'EMD')
    emdComps=emd(signal,'Display',0);%run EMD
    %...and take the first component
    centeredSig=emdComps(:,1);
else
    [~,~,meanEnv]=get_envelopes(signal); %extract mean envelope
    centeredSig=signal-meanEnv'; %center the signal around zero
end
normSig=normalize_env_peaks(centeredSig,[],[],5);%Huang normalizatin algo
% remove potentially divergent behaviour at the edges
normSig=edgeCorrect(normSig, 1.2);
PHI=angle(hilbert(normSig,tmpLen));%hilbert phase
PHI=PHI(1:origLen);

PHIc=correct_PHI(PHI,SG_coef,handleNegPHI);%corrected phase

