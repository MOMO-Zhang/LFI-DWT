function [P] = scramble_M(I,X_M)
%��ͼ��I��������
%   I����Ҫ���ҵ�ͼ�����X_M�����Һ��ͼ�����������У�P�����²�����ͼ��
%Ԥ�����ڴ�
[M,N]=size(I); 
P = zeros(M,N);
for flag_M = 1:M
    P(flag_M,:) = I(X_M(flag_M),:);
end
end