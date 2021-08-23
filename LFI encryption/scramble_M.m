function [P] = scramble_M(I,X_M)
%对图像I进行置乱
%   I是需要置乱的图像矩阵，X_M是置乱后的图像行序数排列，P是最新产生的图像
%预分配内存
[M,N]=size(I); 
P = zeros(M,N);
for flag_M = 1:M
    P(flag_M,:) = I(X_M(flag_M),:);
end
end