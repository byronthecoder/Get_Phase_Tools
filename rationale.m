pathIn='C:\Users\lancia\Documents\MATLAB\myTutos\signals\EMA\cl_pata_2.mat';
dataIn=load(pathIn); %load data 
sigIn=dataIn.data(:,20);%get jawY sensor (cf. strmatch('jawY',dataIn.descriptor,'exact') )
sr=dataIn.samplerate;% sampling rate

%low pass filter parameters
hiCut=25; %Max allowed frequency
filtLen=1000;%filter order

%Savitzky-Golay differentiator parameters:
m = 5;      % number of filtered points (FIR order)
M = (m+1)/2;      % output point , M = (m+1)/2; % middle point for odd m
n = 3;      % approximation polynomial order
h = 1/sr; %sample period
SG_coef = SG_calc(m,n,h); %compute coeffs for the sgolay differentiator
handleNegPHI='interp'; % set strategy to adopt in intervals with negative frequency

%%here we lowpass filter the signal

hiCutNorm=hiCut./( sr/2 ); %normalized cut off for the low-pass filter
b1=fir1(filtLen,hiCutNorm,'low'); %filter coeffs
filtSig=filter(b1,1,sigIn); %apply filter to input sig.
filtSig=zscore(filtSig);%at least now it oscillates around 0 like a real sinusoid.
sigInf=filtSig(filtLen/2+1:end); %correct filter delay

[PHI,PHIc]=getPHI(sigInf,sr,m,n,'',handleNegPHI);

%a successful application of Huang et al 2009 normalization requires a signal centered around zero. 
[maxEnv,minEnv,meanEnv]=get_envelopes(sigInf); %extract maximum minimum and mean envelope
centeredSig=sigInf-meanEnv'; %center the signal around zero

figure;
o(1)=subplot(211);
plot(sigInf)
hold on,
plot(maxEnv)
plot(minEnv)
plot(meanEnv)
title('filtered signal and envelopes')
o(2)=subplot(212);

plot(centeredSig)
title('signal minus mean envelope')

%alternatively we can run EMD ...
emdComps=emd(sigInf,'Display',0);
%...and take the first component
emdComp=emdComps(:,1);

% or ceemd by Colominas et al. 2014 (toolbox required)
% addpath([userpath,'\TBXS\package_emd\EMDs']);
% [modes,~]=ceemd_new(sigInf,0.52,500,2000,1);
% emdComp=modes(1,:)';

origLen=length(sigInf); %original length of the signal
tmpLen=2^nextpow2(origLen);%set optimal signal length to go through the Hilbert transform

normSig=normalize_env_peaks(centeredSig,[],[],5);%Huang normalizatin algo
normSig0=normalize_env_peaks(sigInf,[],[],5);%normalize peaks also for the signal which has not been centered around zero for comparison  
normSig1=normalize_env_peaks(emdComp,[],[],5);%normalize peaks also for the first component of the filtered signal's EMD

% remove potentially divergent behaviour at the edges
normSig=edgeCorrect(normSig, 1.2);
normSig0=edgeCorrect(normSig0, 1.2);
normSig1=edgeCorrect(normSig1, 1.2);

%hilbert phase
PHI0=angle(hilbert(normSig0,tmpLen));
PHI0=PHI0(1:origLen);
PHI=angle(hilbert(normSig,tmpLen));
PHI=PHI(1:origLen);
PHI1=angle(hilbert(zscore(normSig1),tmpLen));
PHI1=PHI1(1:origLen);
PHI2=angle(hilbert(zscore(sigInf),tmpLen));
PHI2=PHI2(1:origLen);

%phase correction 
PHIC=correct_PHI(PHI,SG_coef,handleNegPHI);
PHI0C=correct_PHI(PHI0,SG_coef,handleNegPHI);
PHI1C=correct_PHI(PHI1,SG_coef,handleNegPHI);
PHI2C=correct_PHI(PHI2,SG_coef,handleNegPHI);

figure
a(1)=subplot(5,1,1);
plot(sigInf)
title('filtered signal')

a(2)=subplot(6,1,2);
plot(normSig)
line([1,length(sigInf)],[0,0])
hold on;
plot(normSig0);
title('signal normalized with without centering')
a(3)=subplot(6,1,3);
plot(PHI)
hold on
plot(PHIC)
title('phase from centered and normalized sig')
a(4)=subplot(6,1,4);
plot(PHI0)
hold on
plot(PHI0C)
title('phase from normalized (non centered) sig')
a(5)=subplot(6,1,5);
plot(PHI1)
hold on
plot(PHI1C)
title('phase from normalized first component of signal EMD')
a(6)=subplot(6,1,6);
plot(PHI2)
hold on
plot(PHI2C)
title('phase from somoothed sig')
linkaxes(a,'x')
% PHI=wrapTo2Pi(PHI);
   
    