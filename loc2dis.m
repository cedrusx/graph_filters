function dis = loc2dis(loc,adj)
% Generate distant matrix
% loc: a N-by-2 matrix, of which each row contains the x- and y-location 
%   of one node
% adj: N-by-N adjacent matrix of the nodes
% dis: N-by-N distant matrix of the nodes, of which the element (i,j) is
%   0: if i = j
%   d > 0: where d is the true distance between node i and node j
%   Inf: if node i and node j are not adjacent
% 2011-12-29 by snowflurry

N = size(loc,1);
if exist('adj','var')
    dis = (1 - diag(ones(N,1))) ./ (adj + diag(ones(N,1)));
else
    dis = (1 - diag(ones(N,1)));
end
for i = 1:N
    for j = find(dis(i,1:(i - 1)) == 1)
        dis(i,j) = norm(loc(i,:) - loc(j,:));
    end
end
dis = tril(dis);
dis = dis + dis';
        