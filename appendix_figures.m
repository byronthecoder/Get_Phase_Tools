clear all
warning 'off'

addpath('./gen_utils');
addpath(genpath('./get_phase')); 
if exist('./figures1','dir')<1
    mkdir('./figures1')
end

m = 22;      % number of filtered points (FIR order)
ampCoeff=2;
M = (m+1)/2;      % output point , M = (m+1)/2; % middle point for odd m
n = 5;      % approximation polynomial order
threshs=[NaN,NaN];
phaseMethod={'h','cl'};
phaseMethod1={'h','qs'};

nMasks=22;

%%
sr=1000;
dur=1.5;
t=[1:dur*sr]./sr;
omg=12;
phaseY=t.*omg.*2*pi;
y=sin(phaseY);
phaseY1=t*(omg/1.75)*2*pi;
y1=0.8.*sin(phaseY1);
phaseY2=t*(omg/32)*2*pi;
y2=2.*sin(phaseY2);
phaseYOO=t.*3*omg.*2*pi.*(square(t*(omg/6)*2*pi)>0);
yOnOff=sin(t.*3*omg.*2*pi).*(square(t*(omg/6)*2*pi)>0);

YY=y+y1+y2+0.2*yOnOff;

Decomp=emd(YY)';
normSig=demodulateAmp(Decomp(:,1));%Huang normalizatin algo
normSig1=demodulateAmp(Decomp(:,2));%Huang normalizatin algo

[PHI,siftedSig,PHI0,ceteredSig,mask]=getPHImask(YY,sr,m,n,nMasks,[],{'h','h'},threshs);
[PHI1,siftedSig,PHI0,ceteredSig,mask]=getPHImask(YY,sr,m,n,nMasks,[],phaseMethod,threshs);

figure;
[hb, pos] = tight_subplot(2,1,[0.04,0.05],0.1);
axes(hb(1));
plot(YY,'k')
set(gca,'xtick',[])
title('\textbf{Input signal} (AU)','Interpreter','latex')
axes(hb(2));
plot(PHI,'k')
hold on
plot(PHI1,'k:', 'linewidth',1)

set(gca,'xtick',[])
setPhaseAx(hb(2))
title('\textbf{Phase estimate} (rad)','Interpreter', 'latex')

xlabel('Time steps')
print('./figures1/fig_16','-dtiff', '-r600')

%%
artDimN=[16,18,20];%select relevant channels: ttip 16, llip 18, jawY 20 
cycle1=[3150,3350];%boundaries of slow cycle example
cycle2=[2479,2645];%boundaries of fast cycle example
pathIn='./data./cl_pata_2.mat';
dataIn=load(pathIn); %load data 
sigIn=dataIn.data(1:end,:);%(:,artDimN);%get jawY 20 ,ulip 17, llip 18, ttip 16 sensor (cf. strmatch('lev-sY ',dataIn.descriptor,'exact'),'lev-sY','lev-iY ','lang-aY' )
sr=dataIn.samplerate;% sampling rate
sigT=[1:length(sigIn)]./sr;% time steps

%low pass filter parameters
hiCut=20; %Max allowed frequency
filtLen=1000;%filter order
hiCutNorm=hiCut./( sr/2 ); %normalized cut off for the low-pass filter
b1=fir1(filtLen,hiCutNorm,'low'); %filter coeffs
filtSig=filter(b1,1,sigIn); %apply filter to input sig.
filtSig=zscore(filtSig);%at least now it oscillates around 0 like a real sinusoid.
sigInf=filtSig(filtLen/2+1:end,:); %correct filter delay
sigTf=sigT(1:length(sigInf)); %remove exceding samples' time points

sigInf=sigInf(1:end,:);% remove artefacts at the onset
sigTf=sigTf(1:end);% remove artefacts at the onset
myData=[sigInf(:,artDimN(1)),sigInf(:,artDimN(2)),sigInf(:,artDimN(3))];%extract relevant channels

txtY=[-1.2,-2.2,-1.5];
txtX=8;
artLabels={'TTIP','LLIP','JAW'};
fSize=8;

myPhases=zeros(size(myData)); %initialize matrix containing phase values
myPhases1=zeros(size(myData)); %initialize matrix containing phase values

for dataN=1:size(myData,2)
    [PHI,newIMF,PHI0,centeredSig,mask]=getPHImask(myData(:,dataN),sr,m,n,nMasks,ampCoeff,'h',threshs);
    myPhases(:,dataN)=PHI;
    [PHI,newIMF,PHI0,centeredSig,mask]=getPHImask(myData(:,dataN),sr,m,n,nMasks,ampCoeff,phaseMethod,threshs);
    myPhases1(:,dataN)=PHI;
end

figure;
set(gcf, 'position',[83 198 860 294])
[hb, pos] = tight_subplot(4,1,[0.04,0.05],0.1);
axes(hb(1))
plot(sigTf,myData(:,3),'k')
ylabel('AU')
set(gca, 'xtick',[])
xlabel('')
for i=1:3
    axes(hb(i+1))
    plot(sigTf,myPhases(:,i), 'k');
    hold on
    plot(sigTf,myPhases1(:,i),'k:', 'linewidth',1);
    ylabel('rad')
    if i<3
        set(gca,'xtick',[]);
        xlabel('');
    else
        xlabel('Time (sec)');
    end
end
% myPhases(:,1)
% h=stackedplot(sigTf,[myData(:,3),myPhases(:,:)],'DisplayLabels',{'cm' 'rad' 'rad' 'rad' } ,'color','k');
% ax = findobj(h.NodeChildren, 'Type','Axes');
% set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')

print('./figures1/fig_17','-dtiff', '-r600')

%%

        sr=5000;% set new signal sr 
        pathIn='./data/eggSig.wav';
        [audioIn,oldSr]=audioread(pathIn);% load data
        mySig=resample(audioIn(:,2),sr,oldSr);% resample at new sr
        mySigT=[1:length(mySig)]./sr; % build time samples vector
        
        [PHI,newIMF,PHI0,centSig,mask]=getPHImask(mySig,5000,m,n,nMasks,ampCoeff,{'h','h'},threshs);
        [PHI1,newIMF,PHI0,centSig,mask]=getPHImask(mySig,5000,m,n,nMasks,ampCoeff,phaseMethod,threshs);       
        figure;
        set(gcf, 'position', [100 255 1000 362])

        [aa,aapos]=tight_subplot(3,1,0.08,0.08); %initialize subplots
        axes(aa(1));
        plot(mySigT,zscore(mySig),'k'); %plot input signal
        set(gca,'xtick',[])
        title('\textbf{Electroglottographic signal} (ZSU)','interpreter','latex','fontsize',12)
        axes(aa(2))
        plot(mySigT,demodulateAmp(newIMF), 'k');
        title('\textbf{Demodulated main oscillatory component}','interpreter','latex','fontsize',12)
        set(gca,'xtick',[])
        axes(aa(3))
        plot(mySigT,PHI,'k');% plot phase computed via CL
        hold on;
        plot(mySigT,PHI1,'k:', 'linewidth', 1)% plot Hilbert phase 

        title('\textbf{Phase} (rad)','interpreter','latex','fontsize',12)
        linkaxes(aa,'x')
        xlim([ 0.2,0.3])
        xlabel('Time (sec)')
        print('./figures1/fig_18','-dtiff', '-r600')