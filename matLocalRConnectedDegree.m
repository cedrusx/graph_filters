function [A, Loc, rmin, rmax] = matLocalRConnectedDegree(N,degree)
% 生成随机局部连接网络的连接矩阵
% N个节点在1*1空间内均匀分布
% 生成的网络平均度为degree(若N*degree为奇数则略小于degree)
% 生成的网络最大通信半径属于区间(rmin, rmax)
% 确认生成网络是联通的，否则重新生成
% 需要Bioinformatics Toolbox中的graphshortestpath函数
% (toolbox\bioinfo\graphtheory)

while 1
    % global IX IY;
    Loc = random('unif',0,1,N,2); 
    A = zeros(N, N);
    dis = loc2dis(Loc, ones(N, N));
    dis = dis + tril(ones(N, N)) * 2;
    for n = 1:(degree * N / 2)
        [i, j] = find(dis == min(min(dis)));
        i = i(1); j = j(1);
        A(i, j) = 1;
        A(j, i) = 1;
        dis(i, j) = Inf;
    end
    if min(sum(A)) > 0 && ...
            ~isinf(max(graphshortestpath(sparse(A),1,'directed',false)))
        % make sure there is no isolated island in the network
        rmin = dis(j, i) - 2;
        rmax = min(min(dis));
        return;
    end
end