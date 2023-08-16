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

winLen=0.5;
winStep=winLen/2;
wLenFr=ceil(winLen.*sr);
wStepFr=ceil(winStep.*sr);

%%here we lowpass filter the signal

hiCutNorm=hiCut./( sr/2 ); %normalized cut off for the low-pass filter
b1=fir1(filtLen,hiCutNorm,'low'); %filter coeffs
filtSig=filter(b1,1,sigIn); %apply filter to input sig.
filtSig=zscore(filtSig);%at least now it oscillates around 0 like a real sinusoid.
sigInf=filtSig(filtLen/2+1:end); %correct filter delay

[PHI,PHIc]=getPHI(sigInf,sr,m,n,'EMD',handleNegPHI);

nPts=(length(sigIn)-(wLenFr/2))./wStepFr;
Ts=winLen/2+[1:nPts].*winStep;
newPHI=zeros(length(PHI),1);
allFreqs=zeros(length(Ts),1);
for i=1:nPts
    startPt=(i-1)*wStepFr+1;
    endPt=startPt+wLenFr-1;
    endPt=min([endPt,length(sigInf)]);
    if endPt<=startPt
        continue
    end
    thisChunk=sigInf(startPt:endPt);
    thisPhChunk=PHI(startPt:endPt);
    thisPHIu=unwrap(thisPhChunk);
    %get average frequency
    allFreqs(i)=(thisPHIu(end)-thisPHIu(1))/(2*pi*length(thisChunk));
end
t=[1:length(PHI)]./sr;
timeVarFreq=interp1(Ts,allFreqs,t,'pchip','extrap');
mask=sin(wrapToPi(2*pi*cumsum(timeVarFreq)));


[IMF,maskOmega]=mask_sift(sigInf,mask,[],24);%get signal via masked sifting
IMF=normalize_env_peaks(IMF-mean(IMF),[],[],5);%Huang normalizatin algo
origLen=length(IMF);
nextLen=2^nextpow2(origLen);
newPHI=angle(hilbert(IMF));
tmpPHI=tmpPHI(1:origLen);


PHIc=correct_PHI(tmpPHI,SG_coef,handleNegPHI);%corrected phase

figure
a(1)=subplot(211)
plot(sigInf)
a(2)=subplot(212)
plot(newPHI)
hold on;
plot(PHIc)

linkaxes(a,'x')