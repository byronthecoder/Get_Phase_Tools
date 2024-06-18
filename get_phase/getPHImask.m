function [newPHI,IMF1,PHI,centeredSig,mask ]=getPHImask(signal,sr,m,n,nMasks,ampCoeff,quadMethod,threshs,sndCorrSp)
% input: 
% signal: input signal
% sr: (positive integer) sampling rate
% m: (positive integer; optional, default: 16) number of filtered points for Savitzky-Golay differentiator 
% n: (positive integer; optional, default: 5) polinomial order of the differentiator
% nMasks: (positive integer; optional, default 22)number of masks used to extract the independent mode function via masked EMD
% ampCoeff:(positive real; optional, default: 2) coefficient determining the amplitude of the masks as a
%            proportion of 4*std(signal), which is meant to be a rough estimate of the signal's 
%           range from the std if the signal's values are normally distributed.
% quadMethod: (a string or a cell of two strings, default: {'h','h'}). Method to be use in the 
%            computation of the quadrature signal and phase 'h' stands for Hilbert, 'q' for direct quadrature, 
%            'qs' for direct quadrature interpolated around zeroes and 'cl' for curve length. 
%            If two strings are provided a different method will be adopted in in the first 
%            or the second part of the algorithm. When 'qs' is used the span of the interpolation 
%            interval is by default one on each side of the zero crossing,
%            unless differently defined via the parameter 'sndCorrSp'.
% threshs:(scalar or vector of two positive real values close to zero, default: [1E-10, 1E-10]) threshold for refined 
%          amplitude normalization. If to values, different thresholds will be used in the two parts 
%          of the algorithm.
% sndCorrSp:(positive integer: default1) length of the intervals on each side
%           of the quadrature values close to zero that are interpolated in the coputation of the quadrature 
%           signal when using Sandoval and De Leon's (2017) correction.
% output:
% newPHI: instantaneous phase
% IMF1: output of masked sifting
% PHI: initial phase estimte
% centredSig: centred signal,
% mask: mask signal with inital phase set to 0

% Leonardo Lancia 14/10/2023
% mailto: leonardo.lancia@cnrs.fr

if diff(size(signal)) >0% be sure that the signal is a column vector
    signal=signal';
end
if nargin<3 || isempty(m)
    m=16;% number of filtered points (FIR order; keep it even if you do not want warnings).
end
if nargin<4 || isempty(n)
    n=5;% approximation polynomial order
end
if nargin<5 || isempty(nMasks)
    nMasks=22;
end
if nargin<6 || isempty(ampCoeff)
    ampCoeff=[];
end
if nargin<7 || isempty(quadMethod)
    quadMethod={'h','h'};
elseif iscell(quadMethod) && length(quadMethod)==1
    quadMethod={quadMethod{1},quadMethod{1}};
elseif ischar(quadMethod)
    quadMethod={quadMethod,quadMethod};
end
if nargin<8 || isempty(threshs)
    threshs=[NaN,NaN];
elseif length(threshs)==1
    threshs=[threshs,threshs];
end
if nargin <9 || isempty(sndCorrSp)
    sndCorrSp=0;
end
origLen=length(signal); %original length of the signal

M = (m+1)/2;      % output point , M = (m+1)/2; % middle point for odd m
h=1/sr;
g = SG_calc(m,n, 1)'; %compute coeffs for the sgolay differentiator
p=1;

% centring (sifting with one iteration)
options.MAXITERATIONS=1;
options.MAXMODES=1;
options.FIX=1;
[centeredSig,~]=emd(signal,options);
centeredSig=centeredSig(1,:)';

normSig=demodulateAmp(centeredSig,[],threshs(1),5);%refined amplitude normalization

if strcmpi(quadMethod{1},'q')
    PHI=quadAngle(normSig,0); % direct quadrature phase estimation 
elseif strcmpi(quadMethod{1},'qs')
    PHI=quadAngle(normSig,max(1,sndCorrSp)); % direct quadrature phase estimation with Sandoval and de Leon (2017) correction 
elseif strcmpi(quadMethod{1},'h')
    PHI=wrapTo2Pi(unwrap(angle(hilbert(normSig)))); % hilbert  estimation 
elseif strcmpi(quadMethod{2},'cl')
    [ ~,embedding]=quadAngle(normSig,sndCorrSp);
    NV=[1,0];
    shiftAngle=atan2(NV(2),NV(1));%acos(min(1,max(-1, u(:).' * v(:) / norm(u) / norm(v) )));   
    PHI=co_distproto(embedding', NV');% curve length quadrature phase estimation 
    PHI =wrapTo2Pi(PHI-shiftAngle);% correct for initial phase shift due to the chosen Ponicaré surface sect. 
else
    error('quadMethod can only be ''q'', ''h'', ''qs''  or ''cl''')
end

newPHI=PHI;
%estimation of local frequancies
% sigAmp=get_amp(centeredSig);
% sigFreq1=get_omega(newPHI,sr,m,g,sigAmp./range(sigAmp));
sigFreq=get_omega(newPHI,sr,m,g);

t=[1:length(newPHI)]./sr; %time stamps of the phase signal
if ~isempty(ampCoeff)
    ampCoeff=ampCoeff.*4.*std(signal);
end
[IMF1,~,mask]=mask_sift(signal,sigFreq,[],nMasks,ampCoeff);%get signal via masked sifting
IMF=demodulateAmp(IMF1,[],threshs(2),5);%refined amplitude normalizatin algo

if strcmpi(quadMethod{2},'q')
	newPHI=quadAngle(IMF,0); % direct quadrature phase estimation 
elseif strcmpi(quadMethod{2},'h')
    newPHI=wrapTo2Pi(unwrap(angle(hilbert(IMF)))); % hilbert  estimation 
elseif strcmpi(quadMethod{2},'qs')
    newPHI=quadAngle(IMF,max(1,sndCorrSp)); % direct quadrature phase estimation with Sandoval and de Leon (2017) correction 
elseif strcmpi(quadMethod{2},'cl')
    [~,embedding]=quadAngle(IMF,sndCorrSp);
    NV=[1,0];% normal to the Poincare section
    shiftAngle=atan2(NV(1),NV(2));%(argument are inverted to get the direction of the Poincaré section)%acos(min(1,max(-1, u(:).' * v(:) / norm(u) / norm(v) )));   
    newPHI=co_distproto(embedding', NV');% curve length quadrature phase estimation 
    newPHI =wrapTo2Pi(newPHI-shiftAngle);% correct for initial phase shift due to the chosen Ponicaré surface sect. 
else
    error('quadMethod can only be ''q'', ''h'', ''qs'' or ''cl''')
end



