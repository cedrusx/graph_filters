function [A, Loc, rmin, rmax] = matLocalRConnectedDegree(N,degree)
% ��������ֲ�������������Ӿ���
% N���ڵ���1*1�ռ��ھ��ȷֲ�
% ���ɵ�����ƽ����Ϊdegree(��N*degreeΪ��������С��degree)
% ���ɵ��������ͨ�Ű뾶��������(rmin, rmax)
% ȷ��������������ͨ�ģ�������������
% ��ҪBioinformatics Toolbox�е�graphshortestpath����
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