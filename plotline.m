function h = plotline(varargin)
% Plot a set of lines with the same x-data
% e.g. input multiple y-data
% >> figure; plotline(x,y1,y2,y3);
% or input an matrix Y containing each y-data as columns (same as plot())
% >> figure; plotline(x,[y1,y2,y3]); 
% or input an cell array 
% >> figure; plotline(x,{y1,y2,y3});

spec = {'-sr','-^c','-ob','-dk'};
if iscell(varargin{2})
    ys = varargin{2};
elseif min(size(varargin{2})) > 1
    ys = num2cell(varargin{2},1);
else
    ys = varargin(2:end);
end
if length(ys) > length(spec)
    error('I cannot plot so many lines!');
end
h = zeros(length(ys),1);
for k = 1:length(ys)
    h(k) = plot(varargin{1},ys{k},spec{k},'LineWidth',1.5);
    hold on;
end