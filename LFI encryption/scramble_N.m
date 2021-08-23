function [P] = scramble_N(I,X_N)
%对图像I进列置乱
%   I是需要置乱的图像矩阵，X_N是置乱后的图像列序数排列，P是最新产生的图像
%预分配内存
[M,N]=size(I); 
P = zeros(M,N);
for flag_N = 1:N
    P(:,flag_N) = I(:,X_N(flag_N));
end
end
