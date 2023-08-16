function [IMF,maskOmega]=mask_sift(x,maskOmega,maskAmpCoeff,nPhases,frqEst)
 
if nargin<5 || isempty(frqEst)
    frqEst='ZC';
end
if nargin<4 || isempty(nPhases)
    nPhases=4;
end
if nargin<3 || isempty(maskAmpCoeff)
    maskAmpCoeff=std(x);
end
if nargin<2 || isempty(maskOmega) || (length(maskOmega)==1 && maskOmega==0)
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    IMF0=emd(x,'MaxNumIMF',1,'Display',0);
    IMF0=normalize_env_peaks(IMF0);
    if strcmp(frqEst,'phi')
         IPH=angle(hilbert(IMF0));
        Ifreq=diff(unwrap(IPH));
        IPH(Ifreq<0)=NaN;
%         IPH=interp_NAN(IPH,'nearest');
        maskOmega=nanmean(diff(unwrap(IPH))/(2*pi));%normalized freq      
    elseif length(maskOmega)==1
         maskOmega=ceil(length(zci(IMF0))/2)/length(IMF0);%normalized freq
    end
end
phases=linspace(0,2*pi,nPhases+1);
phases=phases(1:end-1);
if length(maskOmega)==length(x)
    masks=maskOmega'+phases;
else
    masks=range(x).*sin([0:length(x)-1].*maskOmega*2*pi+phases')';
end
tmpIMFS=zeros(size(masks));
for i=1:size(masks,2)
    tmpIMFS(:,i)=emd(x+masks(:,i),'MaxNumIMF',1,'Display',0);
end
IMF=sum(tmpIMFS-masks,2)./nPhases;


