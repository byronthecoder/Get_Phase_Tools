function sigIn=normalize_cycle_amp(sigIn,symmetric,threshNorm,maxIterN, symN)
%    normalize envelope peaks iteratively following a modified version of Huang et al. (2009) 
%    (preventing the occurrence of local extrema smaller than 1 and larger than -1) and using pchip
%     interpolation to obtain the envelope connecting the abs. values of local extrema.
%    Boundary conditions of the local maxima sequence are computed following 
%    Rilling 2007 (see function in '.\utils\boundary_conditions.m').
%
%     input:
%         symmetric (default: 1) : if it is 0 we compute the signal minimum and maximum envelopes by 
%         interpolating the local low and high extrema. The mean of these two envelopes 
%         is also computed and it is subtracted by the signal 
%     
%         threshNorm: threshold to determine if peaks are different from 1/-1
%         maxIterN: maximum number of iterations
%         symN: number of extrtema to be added at the beginning and at the
%               end before interpolating.
%     output:
%         normalized envelope
% 
%         sigIn=copy.copy(sig)
% 

% Leonardo Lancia 14/10/2023
% mailto: leonardo.lancia@cnrs.fr

if diff(size(sigIn))>0 % be sure that the input signal is a column vector
    sigIn=sigIn';
end

if nargin <2 || isempty(symmetric)
    symmetric=1;
end
if nargin <3 || isempty(threshNorm) || isnan(threshNorm)
    threshNorm=1e-10;
end
if nargin<4 || isempty(maxIterN)
    maxIterN=5;
end

if nargin<5 || isempty(symN)
    symN=2;
end
if symmetric==0 %if the signal is not centred on 0apply centring (one sifting iteration)
    options.MAXITERATIONS=1;
    options.MAXMODES=1;
    [sigIn,~]=emd(sigIn,options);
    sigIn=sigIn(1,:)';
end


newIdxs=[1:length(sigIn)]'; %input signal's time steps 

iter_val=1; %intitialize error value
niter=0; % current iteration number
while niter<maxIterN && iter_val>threshNorm
    niter=niter+1; %update iteration number
    
    [~,pksLoc]=findpeaks(sigIn); % find local maxima
    [~,valLoc]=findpeaks(-sigIn); %find local minima
    if length(pksLoc)+length(valLoc)<3 %if less than three extrema break and return
        sigIn=sigIn/range(sigIn);
        return 
    end    
    %set boundary conditions via Rilling's function
    [tMin,tMax,sigMin,sigMax] = boundary_conditions(valLoc',pksLoc',newIdxs',sigIn',sigIn',symN);
    %find local extrema order and sort their locations
    [tExtr,extrOrd]=sort([tMin,tMax]);
    allExtr=[sigMin,sigMax];
    allExtr=abs(allExtr(extrOrd));%sort local extrema absolute values according ot their order of appearence

    refSig= interp1(tExtr,allExtr,newIdxs, 'pchip'); %interpolate local extrema absolute values

	iter_val = abs(sum(refSig) - max(size(refSig))) ; % compute mean difference between envelopes values from last two iterations
    sigIn=sigIn./refSig;% normalize signal's amplitude

end