function quadAngle = quadAngle(nIMF)
%
% INPUT:
%       nIMF     - normalized IMFs to calculate analytical signal on
%                  (an 1D data) 
%
% OUTPUT:
%       quadrature  -the analytical signal
%                    A complex array with actual data stored in real part 
%                    and  the phase-shifted equivalent in the imaginary part
%
% NOTE: 
%       Quadrature- a new method to calculate instantaneous frequency propose by Norden E. Huang
%                   There are 3 steps .
%                    1.normalize the IMF, this helps IMF become AM/FM disjointed 
%                    2.form the analytical signal,the normalized-IMF(x), be the real part
%                                                 the quadrature of normalized-IMF(x)- sqrt(1-xx) , be the imaginary part
%                    3.calculate the nstantaneous frequency,use real/imaginary part to find phase angle,and the time derivative is I.F
%       this code is some part of the quadrature calculations
%
%       noh.m------use Normalized IMF to compose the analytical signal using quadrature calculation
%       noh() - creates the quadtrature signal without the Hilbert transform
%
% References:   
%  N. E Huang (2008),NCU hht class lecture 
%  6 Instantaneous Frequency.ppt
%
% code writer: Karin Blank (karin.blank@nasa.gov) for Norden Huang (5/6/05),code originally from nfame.m
% footnote:S.C.Su 2009/04/18
%
% There are only one loop in this code dealing with 1D IMF.
%  1.check input,flip signal if necessary
%  2.Calculate the quadrature value--loop start 
%     3.add mask-for positive/negative signal  
%     4.the calculation of quadrature value   
%     5.multiplying by mask flips data in 3rd & 4th quadrature   
%  2.Calculate the quadrature value--loop end      
%     
%     
% Association: no
% 1. this function ususally used for calculation for I.F
%    the second part-form the analytical signal
%    this code only dealing with forming analytical signal
% 
% Concerned function: epnormal,calcqf
%                     the others are matlab functions.  
% 

%1.check input,flip signal if necessary
%----- Get the dimension
[npt,ncol] = size(nIMF);

%----- Initialize and flip data if needed 
flipped=0;
if (ncol > npt)
    flipped=1;
    nIMF=nIMF';
    [npt,ncol] = size(nIMF);
    
end

%2.Calculate the quadrature value--loop start
 
 
for i=1:ncol;
        data = nIMF(:,i);
       
%3.add mask-for positive/negative signal    
        %creates a mask for the signal 
        %the 1st & 2nd Q = 1 
        %and the 3rd & 4th Q = -1
        mask = ((diff(data)>0) .* -2) + 1;
        mask(end+1) = mask(end);
    
%4.the calculation of quadrature value
        y = real(sqrt(1-data.^2));
    
%5.multiplying by mask flips data in 3rd & 4th quadrature
        q = y .* mask;
        quadrature(:,i) = complex(data, q);
        clear mask y q
end
quadAngle=angle(quadrature);
 if flipped==1
    quadAngle = quadAngle';
    end
%2.Calculate the quadrature value--loop end
%end of noh.m   