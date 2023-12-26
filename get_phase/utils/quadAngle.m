function quadAngle = quadAngle(nIMF,hilbQuad, thresh)
%
% INPUT:
%       nIMF     - matrix containing normalized IMFs to calculate analytical signal on
%       hilbQuad - toggle : 0 (default),1. If hilbQuad==1 then hilbert
%                  quadrature is used (mask = sign(hilbert(data)) otherwise direct  
%                  quadrature (mask= ((diff(data)>0) .* -2) + 1; mask(end+1) = mask(end);)
%       thresh   - positive float <1. If hilbQuad is equal to 0 and thresh>0 elements of 
%                  the quadrature signal whose absolute values are smaller than thresh are
%                  replaced by interpolated values (pchip interpolationis used).
   
% OUTPUT:
%       quadAngle: array of angles obtained from data and quadrature
%                 signals
%

if nargin<2 || isempty(hilbQuad)
    hilbQuad=0;
end
if nargin<3 || isempty(thresh)
    thresh=0;
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
        if hilbQuad==0 %quadrature from signal direction of change
            
            mask = [((diff(data)>0) .* -2) + 1;0];
            mask(end) = mask(end-1);
        else %quadrature from sign of Hilbert transform
            P=hilbert(data);
            mask=sign(imag(P));
        end
        % quadrature value
        y = real(sqrt(1-data.^2));
    
        %flip data in 3rd & 4th quadrants
        q = y .* mask;
        if thresh>0
            q(abs(q)<0.4)=NaN;
            q=interp_NAN(q);
        end
        %store column
        quadrature(:,i) = complex(data, q);
end
quadAngle=wrapTo2Pi(unwrap(angle(quadrature)));%compute angles

%flip data bask if needed 

if flipped==1
    quadAngle = quadAngle';
end
