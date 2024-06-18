function [quadAngle,quadrature] = quadAngle(nIMF,sandLeonCorr,thresh)
%
% INPUT:
%       nIMF     - matrix containing normalized IMFs to calculate analytical signal on
%       hilbQuad - toggle : 0 (default),1: If hilbQuad==1 then to obtain the mask hilbert
%                  quadrature is used (mask = sign(hilbert(data)) otherwise direct  
%                  quadrature (mask= ((diff(data)>0) .* -2) + 1; mask(end+1) = mask(end);)
%       sandLeonCorr - semi positive integer: if >0 Sandoval and De Leon (2017) correction is applied. 
%            points smaller than thresh of the quadrature signal, as well as all points in the 
%            intervals sourronding them and whose length is sandLeonCorr are substituted 
%            by interpolated values.
%        thresh= threshold below which quadrature signals are considered equal to 0
%            for the purposes of the Sandoval and De Leon (2017) correction.
                   
% OUTPUT:
%       quadAngle: array of angles obtained from data and quadrature
%                 signals
%


if nargin<2 || isempty(sandLeonCorr)
    sandLeonCorr=0;
end
if nargin<3 || isempty(thresh)
    thresh=1e-2;
end
if sandLeonCorr>0
    if mod(sandLeonCorr,2)==0
        halfCorrInt=sandLeonCorr/2;
    else
        halfCorrInt=(sandLeonCorr-1)/2;
    end
end

[npt,ncol] = size(nIMF);

%flip data if needed 
flipped=0;
if (ncol > npt)
    flipped=1;
    nIMF=nIMF';
    [npt,ncol] = size(nIMF);
end

%initialize quadrature array
quadrature=zeros(size(nIMF,1),ncol);
%loop over columns
for i=1:ncol
        data = nIMF(:,i);
       %create mask
        
        mask = ((gradient(data)>0) .* -2) + 1;
        mask(mask==0)=-1;

        % quadrature value
        y = real(sqrt(1-data.^2));
    
        %flip data in 3rd & 4th quadrants
        q = y .* mask;
        if sandLeonCorr>0
            AA=find(abs(q)<thresh);%sort(zci(-q));
            AA0=AA-(halfCorrInt);
            AA1=AA+halfCorrInt;

            tmpMat=arrayfun(@(x,y,z) linspace(x,y,y-x+1),AA0,AA1,'UniformOutput',false);
            tmpMat=cell2mat(tmpMat)';
            tmpMat=tmpMat(:);
            tmpMat=tmpMat(tmpMat>0 & tmpMat<=length(q));
            tmpMat=unique(tmpMat);
            q(tmpMat)=nan;
            q=interp_NAN(q);
        end
        
        %store column
        quadrature(:,i) = complex(data, q);
end
quadAngle=wrapTo2Pi(unwrap(angle(quadrature)));%compute angles

%flip data bask if needed 

if flipped==1
    quadAngle = quadAngle';
    quadrature=quadrature';
end
