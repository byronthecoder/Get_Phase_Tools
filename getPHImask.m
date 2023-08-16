function [PHI,PHIc,newPHI,newPHIc]=getPHImask(signal,sr,m,n,winLen,winStep,nMasks, handleNegPHI)

%input: 
% signal: input signal
% sr: sampling rate
% m: number of filtered points for Savitzky-Golay differentiator 
% n: polinomial order of the differentiator
% winLen: duration (in sec) of the window used to esitmate the instantaneous frequency 
% winStep: the hop size for the frequency estimation
% nMasks: number of masks used to extract the independent mode function via masked EMD
% handleNegPHI: parameter handling negative frequency correction. Can be interp or a scalar, in
%   the former case the values are interpolated otherwise they are replaced
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
    winLen=0.5;
end
if nargin<6
    winStep=winLen/2;
end
if nargin<7
    nMasks=8;
end
if nargin<7
    handleNegPHI='interp';
end


M = (m+1)/2;      % output point , M = (m+1)/2; % middle point for odd m
h = 1/sr; %sample period
SG_coef = SG_calc(m,n,h); %compute coeffs for the sgolay differentiator
wLenFr=ceil(winLen.*sr);
wStepFr=ceil(winStep.*sr);
origLen=length(signal); %original length of the signal
tmpLen=2^nextpow2(origLen);%set optimal signal length to go through the Hilbert transform

% 
% emdComps=emd(signal,'Display',0);%run EMD
% %...and take the first component
% centeredSig=emdComps(:,1);

[~,~,meanEnv]=get_envelopes(signal); %extract maximum minimum and mean envelope
centeredSig=signal-meanEnv'; %center the signal around zero


normSig=normalize_env_peaks(centeredSig,[],[],5);%Huang normalizatin algo
% remove potentially divergent behaviour at the edges
normSig=edgeCorrect(normSig, 1.2);
PHI=angle(hilbert(normSig,tmpLen));%hilbert phase
PHI=PHI(1:origLen);

PHIc=correct_PHI(PHI,SG_coef,handleNegPHI);%corrected phase

nPts=(length(signal)-(wLenFr/2))/wStepFr;
Ts=winLen/2+[1:nPts].*winStep;

allFreqs=zeros(length(Ts),1);
for i=1:nPts
    startPt=(i-1)*wStepFr+1;
    endPt=startPt+wLenFr-1;
    endPt=min([endPt,length(signal)]);
    if endPt<=startPt
        continue
    end
    thisChunk=signal(startPt:endPt);
    thisPhChunk=PHIc(startPt:endPt);
    thisPHIu=unwrap(thisPhChunk);
    %get average frequency
    allFreqs(i)=(thisPHIu(end)-thisPHIu(1))/(2*pi*length(thisChunk));
end
t=[1:length(PHI)]./sr;
timeVarFreq=interp1(Ts,allFreqs,t,'pchip','extrap');
mask=sin(wrapToPi(2*pi*cumsum(timeVarFreq)));


[IMF,~]=mask_sift(signal,mask,[],nMasks);%get signal via masked sifting
IMF=normalize_env_peaks(IMF-mean(IMF),[],[],5);%Huang normalizatin algo
IMF=edgeCorrect(IMF, 1.2);
origLen=length(IMF);
nextLen=2^nextpow2(origLen);
newPHI=angle(hilbert(IMF,nextLen));
newPHI=newPHI(1:origLen);
newPHIc=correct_PHI(newPHI,SG_coef,handleNegPHI);%corrected phase
