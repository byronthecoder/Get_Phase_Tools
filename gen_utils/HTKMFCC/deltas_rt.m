function  [d,newFiltStCond] = deltas_rt(x, w, varargin)
% D = deltas(X,W)  Calculate the deltas (derivatives) of a sequence
%    Use a W-point window (W odd, default 9) to calculate deltas using a
%    simple linear slope.  This mirrors the delta calculation performed 
%    in feacalc etc.  Each row of X is filtered separately.
% 2003-06-30 dpwe@ee.columbia.edu


if nargin>1 && ~isempty(varargin)
   filtStCond= varargin{1};
else
    filtStCond=[];
end
    

% Define window shape
nPastVals = w - 1;
hlen=floor(nPastVals/2);
win = hlen:-1:-hlen;%nPastVals:-1:0;
%pad if necessary
% % if size(x,2)<w
% %     x=[repmat(x(:,1),1,nPastVals),x];
% % end
% Apply the delta filter
if filtStCond==0
    [d,newFiltStCond] = filter(win, 1, x, [], 2);  % filter along dim 2 (rows)
else
    [d,newFiltStCond] = filter(win, 1, x, filtStCond, 2); 
end
