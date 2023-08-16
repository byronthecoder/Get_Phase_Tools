function d = distance(a,b,thisnorm,varargin)
%this function is built around the code by R. Bunschoten (see below). if thisnorm =='Euclidean' (default), euclidean thisnorm is computed. 
%Otherwise, if thisnorm=='Maximum', superior thisnorm is computed. 
%the matrices' dimensions chould be on different rows.

% Input: matrix a
%	     matrix b
%          thisnorm ('Euclidean' or 'Maximum')
% Output: Distance matrix

% it follows the comment to the original code for euclidean distance 
% DISTANCE - computes Euclidean distance or maximum thisnorm distance matrix

% E = distance(A,B)
%
%    A - (DxM) matrix 
%    B - (DxN) matrix
%
% Returns:
%    E - (MxN) Euclidean distances between vectors in A and B
%
%
% Description for Euclidean distance: 
%    This fully vectorized (VERY FAST!) m-file computes the 
%    Euclidean distance between two vectors by:
%
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Example : 
%    A = rand(400,100); B = rand(400,200);
%    d = distance(A,B);

% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Oct 29 16:35:48 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3
% Thanx    : Nikos Vlassis

% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.
%
% december 2012: Leonardo Lancia added superior norm

if (nargin < 2)
   error('Not enough input arguments');
end

if (nargin < 3)
   thisnorm='Euclidean';
end

if (size(a,1) ~= size(b,1))
   error('A and B should be of same dimensionality');
end

%LL: if vectors it corrects for errors around 0 due to sampling rate:
if sum(size(a)>1)==1 && sum(size(b)>1)==1
    if diff(size(a))<0
        a=a';
    end
    if diff(size(b))<0
        b=b';
    end
end
%     d=repmat(a,1,length(a))-repmat(b',length(b),1);
%     for i=1:size(d,2)
%         ocross=find((d(1:end-1,i)<0 & d(2:end,i)>0) | (d(1:end-1,i)>0 & d(2:end,i)<0) );
%         pathacross=[ocross-1,ocross,ocross+1];
%         pathacross(pathacross==0)=1;
%         for j=1:length(ocross)
%             [~,goodidx]=min(d(pathacross(j,:),i));
%             d(pathacross(j,goodidx),i)=0;
%         end
%     end
%     d=sqrt(d.^2)';
% else

if ~isempty(varargin) && size(a,2)*size(b,2)>varargin{1}
    blockSize=varargin{1};
    inc=round(min([size(b,2)/4,size(a,2)/4,blockSize]));startA=1;
    lenB=size(b,2);lenA=size(a,2);
    d=zeros(lenB,lenA);
    while startA<=size(a,2)
        endA=min([startA+inc,lenA]);
        startB=1;
        while startB<=size(b,2)
            endB=min([startB+inc,lenB]);
            thisA=a(:,startA:endA);
            thisB=b(:,startB:endB);
            thisDist=distance(thisA,thisB,thisnorm);
            d(startB:endB,startA:endA)=thisDist;
            startB=endB+1;
        end
        startA=endA+1;
    end
    
else

    % LL: switch between the two norms
    if ~isempty(strmatch(thisnorm,'Maximum') ) &  sum(size(a)>1)>1%superior thisnorm
        a=a';b=b';
        px=permute(a,[1,3,2]);
        py=permute(b,[3,1,2]);
        s1 = px(:,ones(1,size(b,1)),:) - py(ones(1,size(a,1)),:,:);
        d = max(abs(s1),[],3);d=d';
    else    %euclidean norm
%         ~isempty(strmatch(thisnorm,'Euclidean')) 
        aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b; 
        d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
        d=d';
    end
end
