
%Savitzky-Golay differentiator parameters:
m = 5;      % number of filtered points (FIR order)
M = (m+1)/2;      % output point , M = (m+1)/2; % middle point for odd m
n = 3;      % approximation polynomial order
h = 1/sr; %sample period
SG_coef = SG_calc(m,n,h); %compute coeffs for the sgolay differentiator
handleNegPHI='interp'; % set strategy to adopt in intervals with negative frequency

winLen=0.5;
winStep=winLen/2;
nMasks=16;
exampleN=2;
switch exampleN
    case 1

        pathIn='C:\Users\lancia\Documents\MATLAB\myTutos\signals\EMA\cl_pata_2.mat';
        dataIn=load(pathIn); %load data 
        sigIn=dataIn.data(:,20);%get jawY sensor (cf. strmatch('jawY',dataIn.descriptor,'exact') )
        sr=dataIn.samplerate;% sampling rate

        %low pass filter parameters
        hiCut=25; %Max allowed frequency
        filtLen=1000;%filter order
        % %%here we lowpass filter the signal

        hiCutNorm=hiCut./( sr/2 ); %normalized cut off for the low-pass filter
        b1=fir1(filtLen,hiCutNorm,'low'); %filter coeffs
        filtSig=filter(b1,1,sigIn); %apply filter to input sig.
        filtSig=zscore(filtSig);%at least now it oscillates around 0 like a real sinusoid.
        sigInf=filtSig(filtLen/2+1:end); %correct filter delay
        [PHI,PHIc,newPHI,newPHIc]=getPHImask(sigInf,sr,m,n,winLen,winStep,nMasks,handleNegPHI);


        figure
        a(1)=subplot(311);
        plot(sigInf)
        a(2)=subplot(312);
        plot(PHI)
        hold on;
        plot(PHIc)
        a(3)=subplot(313);
        plot(newPHI)
        hold on;
        plot(newPHIc)
        linkaxes(a,'x')

    case 2
        srAudio=16000;
        srAmp=1000;
        fileName='FR_F_AB4_2016_10_24_ModuleDiadoco_bababa';
        inTxtFile='textgrid_info_FR_GE.csv';
        spGroup=1;
        chunksXfile=1;
        [AMsigs,AMt,srAmp]=getAmpSigs(fileName,inTxtFile,spGroup,chunksXfile);

        [PHI,PHIc,newPHI,newPHIc]=getPHImask(AMsigs(:,1),srAmp,m,n,winLen,winStep,nMasks,handleNegPHI);

        figure
        b(1)=subplot(211);
        plot(AMsigs(:,1))
        hold on;
        plot(AMsigs(:,2))
        b(2)=subplot(212);
        plot(PHI)
        hold on;
        plot(newPHIc)
        linkaxes(b,'x')

    case 3
        srAudio=16000;
        srAmp=1000;

        pathIn='C:\Users\lancia\Desktop\reading_corpora\output\ivie\c-text-f.wav';
        [signal,sr]=audioread(pathIn);
        signal1=resample(signal,srAudio,sr);
        sr=srAudio;
        segRate=16;
        syllRate=5.5;
        stressRate=3;
        winLenBPass=1000;
        nshift=floor(winLenBPass/2);
        AMsig=abs(hilbert(zscore(signal)));

        lpfilt =fir1(500,[1/(FSsig/2) (0.75*segRate)/(FSsig/2)]);

            AMsig=filter(lpfilt,1,AMsig);
            AMsig=sample_advance(AMsig, floor(500/2), 1e-7);

            % newCentFreq = findOptOScFreq2(AMsig,pars.FSamp,[],[pars.startT,pars.endT]);

            AMsig=resample(AMsig,FSamp,FSsig);


            tAmp=[1:length(AMsig)]./FSamp;

            CF0 = [stressRate/2;stressRate;syllRate;segRate;max(segRate,syllRate*3)];
            edges=CF0(1:end-1)+diff(CF0)/2;

            [CF, bpfs] = MFB_coeffs(edges,FSamp,0);


            AMsigs=zeros(size(AMsig,1),3);
             for n = 1:length(CF)
                l_fil = bpfs(1,n);
                fil_n = bpfs(2:1+l_fil,n);
                nshift = floor(l_fil/2);
                AMfil = filter(fil_n, 1, AMsig);
                AMfil = sample_advance(AMfil, nshift, 1e-7);
                AMsigs(:,n) =AMfil;
             end

        [PHI,PHIc,newPHI,newPHIc]=getPHImask(AMsigs(:,2),srAmp,m,n,winLen,winStep,nMasks,handleNegPHI);

        figure
        b(1)=subplot(211);
        plot(AMsigs(:,2))
        b(2)=subplot(212);
        plot(PHI)
        hold on;
        plot(newPHIc)
        linkaxes(b,'x')
end