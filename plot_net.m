function plot_net(loc, adj)
% plot the figure of a network
% loc: a matrix contains the locations of nodes, with x- and y-location
%   of one node in each row
% adj: adjacent matrix of the nodes

hold on;
N = size(loc,1);
[is,js] = find(tril(adj));
for k = 1:length(is)
    i = is(k); j = js(k);
    plot([loc(i,1),loc(j,1)],[loc(i,2),loc(j,2)],'-m');
end
plot(loc(:,1),loc(:,2),'.b');
hold off;