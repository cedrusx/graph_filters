function [A, Location] = matLocalRConnected(N,d)
% ��������ֲ�������������Ӿ���
% N���ڵ���1*1�ռ��ھ��ȷֲ�
% ���ڵ�䵱ŷ�Ͼ��벻����dʱ���ڱ�����
% ȷ��������������ͨ�ģ�������������
% ��ҪBioinformatics Toolbox�е�graphshortestpath����
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