function [P] = scramble_N(I,X_N)
%��ͼ��I��������
%   I����Ҫ���ҵ�ͼ�����X_N�����Һ��ͼ�����������У�P�����²�����ͼ��
%Ԥ�����ڴ�
[M,N]=size(I); 
P = zeros(M,N);
for flag_N = 1:N
    P(:,flag_N) = I(:,X_N(flag_N));
end
end
