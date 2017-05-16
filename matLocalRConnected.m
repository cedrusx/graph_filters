function [A, Location] = matLocalRConnected(N,d)
% 生成随机局部连接网络的连接矩阵
% N个节点在1*1空间内均匀分布
% 两节点间当欧氏距离不大于d时存在边连接
% 确认生成网络是联通的，否则重新生成
% 需要Bioinformatics Toolbox中的graphshortestpath函数
% (toolbox\bioinfo\graphtheory)

while 1
    % global IX IY;
    IX=random('unif',0,1,N,1); 
    IY=random('unif',0,1,N,1);
    dd=d*d;
    A=zeros(N,N);
    for i=1:N
        for j=(i+1):N
            h=IX(i)-IX(j);
            l=IY(i)-IY(j);
            if h^2+l^2<=dd
                A(i,j)=1; A(j,i)=1;
            end
        end
    end
    if min(sum(A)) > 0 && ~isinf(max(graphshortestpath(sparse(A),1,'directed',false)))
        % make sure there is no isolated island in the network
        Location = [IX IY];
        return;
    end
end