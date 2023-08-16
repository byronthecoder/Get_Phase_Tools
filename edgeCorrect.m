function sigIn=edgeCorrect(sigIn, thresh)
% remove potentially divergent behaviour at the edges by setting to nan and
% then interpolating values whose magnitude exceeds a threshold
% 
% input: 
%     sigIn: signal to be processed
%     thresh: threshold for absolute values
% output: corrected signal
% 
sigIn(abs(sigIn)>thresh)=nan;
if sum(isnan(sigIn))>0
    myIdxs=find( ~isnan(sigIn));
    myVals=sigIn(myIdxs);
    [sigIntmp,sigIntmplocs]=SetBoundCond(myVals,myIdxs,length(sigIn));
    sigIn=interp1(sigIntmplocs,sigIntmp,[1:length(sigIn)], 'pchip');
end
