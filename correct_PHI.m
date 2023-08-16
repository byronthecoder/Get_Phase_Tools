function PHI=correct_PHI(PHI,SG_coef,handleNegPHI)
%in order to correct potentially artefactual phase values, first we
%estimate the phase difference, then we indentify intervals where the phase
%difference is negative. The endpoint of each identified interval 
% is then moved to the first time step where the phase is higer than the initial
%phase of the interval. values contained in the intervals are set to nan
%and interpolated or are replaced by a constant value set by handleNegPHI
% input: 
% PHI: phase signal
% SG: differentiator coeffs
%handleNegPHI: parameter handling negative frequency correction. Can be interp or a scalar, in
%   the former case the values are interpolated otherwise they are replaced
%   with the scalar

if nargin<3
    handleNegPHI='interp';
else
    if ~ischar(handleNegPHI)
        fillVal=handleNegPHI;
        handleNegPHI='fill';
    end
end
if diff(size(PHI))>0
    PHI=PHI';
end
m=size(SG_coef,2);
PHIu=unwrap(PHI);
badGuys=zeros(length(PHI),1);

dy_est  = filter(SG_coef(2,:),1,PHIu);
ind=m:(length(PHI)-m);
deltaPHI=dy_est(ind);
%     deltaPHI=sgolayfilt(diff(sgolayfilt(unwrap(PHI), 3,5)), 3,5);
badGuys(deltaPHI<0)=1;

badGuys=[0;badGuys];
startBad=find(diff(badGuys)==1)+1;
endBad0=find(diff(badGuys)==-1);

if sum(badGuys)>0
    endBad=zeros(length(startBad),1);
    allPos=[1:length(PHI)]';
    badGuys=zeros(length(PHI),1);
    for nn =1:length(startBad)
        allCandEnd=find( PHI>PHI(startBad(nn)) & allPos>startBad(nn) );
        if ~isempty(allCandEnd)
            endBad(nn)=allCandEnd(1);
        else
            if nn==length(startBad)-1
                endBad(nn)=length(PHI) ;
            else
                endBad(nn)=endBad0(nn)+1;
            end
        end
        badGuys(startBad(nn):endBad(nn))=1;
    end

    if strcmp(handleNegPHI,'interp')
    %         PHI(badGuys==1)=nan;
        PHIu(badGuys==1)=nan;
    %         PHI1=interp_NAN(PHI,'linear');
        PHI=wrapToPi(interp_NAN(PHIu,'linear'));
    %         if sum(diff(PHI1)<0) < sum(diff(PHI2)<0)
    %             PHI=PHI1;
    %         else
    %             PHI=PHI2;
    %         end
    elseif strcmp(handleNegPHI,'fill')   
        PHI(badGuys==1)=fillVal;
    end
end