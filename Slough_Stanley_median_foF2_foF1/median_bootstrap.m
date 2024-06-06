function [median_data, SE_median] = median_bootstrap(data, min_count)

% A function to estimate the standard error in the median of a distribution
% using a bootstrap method - as described in Computational Statistics using
% Matlab.
% 
% Inputs - the data set
% Output - the median and the standard error in the median, SE_median
%
% CJD April 2013
%
% Updated to account for NaNs in the data CJS 2023

if exist('min_count')==0
    min_count = 1;
end

good = find(isnan(data) == 0);
data = data(good);
%data = data(~isnan(data));

B = 100;

if length(data) >= min_count
   median_data = median(data);
   l_data = length(data);
   pseudo_median = zeros(1,B);

   R=randi(l_data,[B,l_data]);
   pseudo_median = median(data(R),2);

% for i=1:B
%     R=randi(l_data, [1,l_data]);
%     pseudo_median(i) = median(data(R));
% end

   mean_pseudo = sum(pseudo_median)/B;

   SE_median = sqrt(( (1/(B-1))*(sum((pseudo_median - mean_pseudo).^2)) ));
else
    median_data = NaN;
    SE_median = NaN;
end

end