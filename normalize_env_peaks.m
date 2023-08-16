function normSig=normalize_env_peaks(sig,symmetric,threshNorm,maxIterN,pksLoc)
%    normalize envelope peaks iteratively following Huang et al. (2009) and using pchip
%     interpolation to obtain the envelope connecting the abs. values of local extrema 
%     input:
%         symmetric (default: 1) : if it is 0 we compute the signal minimum and maximum envelopes by 
%         interpolating the local low and high extrema. The mean of these two envelopes 
%         is also computed and it is subtracted by the signal 
%     
%         threshNorm: threshold to determine if peaks are different from 1/-1
%         maxIterN: maximum number of iterations
%     output:
%         normalized envelope
% 
%         normSig=copy.copy(sig)
% 
normSig=sig;

if nargin <2 || isempty(symmetric)
    symmetric=1;
end
if nargin <3 || isempty(threshNorm)
    threshNorm=1e-10;
end
if nargin<4 || isempty(maxIterN)
    maxIterN=3;
end
if nargin<5
    pksLoc=[];
end
if symmetric==0
    [maxEnv, minEnv,meanEnv]=get_envelopes(normSig);
    normSig=normSig-meanEnv';
end
if isempty(pksLoc)
    [~,pksLoc,proms]=findpeaks(abs(normSig),'MinPeakProminence',0);
%     [~,pksLoc,proms]=findpeaks(abs(normSig),'MinPeakProminence',median(proms)/100);
    extrGiven=0;
else
    extrGiven=1;
end        
    
pksVal=abs(normSig(pksLoc));
if isempty(pksLoc)
    normSig=normSig/range(normSig);
    return 
end

[pksVal,pksLoc]=SetBoundCond(pksVal,pksLoc,length(normSig));
newIdxs=[1:length(normSig)]';
refSig= interp1(pksLoc,pksVal,newIdxs, 'pchip');

iter_val=1;
niter=0;
while niter<maxIterN && iter_val>threshNorm
    niter=niter+1;
    normSig=normSig./refSig;
    [~,pksLoc1,proms]=findpeaks(abs(normSig),'MinPeakProminence',0);
%     [~,pksLoc1,~]=findpeaks(abs(normSig),'MinPeakProminence',median(proms)/100);
    if length(pksLoc1)<2
        normSig=normSig/range(normSig);
        return
    end
    if extrGiven==1
        [pksLoc,~]=get_closest_Vals_Idxs(pksLoc,pksLoc1);

        pksLoc=unique(pksLoc);
    else
        pksLoc=unique(pksLoc1);

    end
    pksVal=abs(normSig(pksLoc));
       
    if length(pksLoc)<2
        return 
    end   
    [pksVal,pksLoc]=SetBoundCond(pksVal,pksLoc,length(normSig));
    refSig= interp1(pksLoc,pksVal,newIdxs, 'pchip');


    iter_val = abs(sum(refSig) - size(refSig,1)) ;
    
end

normSig=normSig./ refSig ;
