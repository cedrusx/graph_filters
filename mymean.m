function [meanx,count] = mymean(x,alpha)
% Revision: 2014.12.7

if nargin < 2
    alpha = Inf;
end
count = zeros(1,size(x,2));
meanx = zeros(1,size(x,2));
for c = 1:size(x,2)
    % Remove NaN
    xc = x(:,c);
    xc = xc(~isnan(xc));
    % Remove outliers
    dis = abs(xc - median(xc)) / mad(xc,1);
    xc = xc(dis <= alpha);
    % Get result
    meanx(c) = mean(xc);
    count(c) = length(xc);
end

% % Remove NaN
% count = sum(~isnan(x));
% y = x;
% y(isnan(x)) = 0;
% m = sum(reshape(y,size(x))) ./ count;
