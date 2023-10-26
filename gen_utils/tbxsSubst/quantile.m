function percentiles=quantileL(data,samplePoints)

%x = 3*randn(1,100)+5;                                                   % Create Data
%h = histogram(x, 20);
%barvals = h.Values;
%edgs = h.BinEdges;
%ctrs = edgs(1:end-1) + (diff(edgs)/2);   

%auc = 100*cumsum(barvals)/sum(barvals)+(0:numel(barvals)-1)*eps*100;    % Cumulative Area
%prctls = [10, 50, 90];                                                  % Desired Percentiles
%prctlctrs = interp1(auc, ctrs, prctls);  

%data = randn(10,3) % 3 sets with 10 samples each
%samplePoints = [0.1, 0.5, 0.9];
percentiles = interp1(linspace(0, 1, size(data,1)), sort(data), samplePoints);  