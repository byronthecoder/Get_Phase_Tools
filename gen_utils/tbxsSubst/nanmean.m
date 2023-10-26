function y = nanmeanL(varargin)
% computes the mean of the first argument once nans have been discarded
y = mean(varargin{:},'omitnan');
