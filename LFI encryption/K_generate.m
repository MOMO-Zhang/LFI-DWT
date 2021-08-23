function [K,Kbin] = K_generate(I)
P_sum=sum(I(:));
P_sum_bin= dec2bin(P_sum);
P_sum_bin = boolean(P_sum_bin-'0');
L_P_sum = length(P_sum_bin);
L_Kbin = 512 - L_P_sum;
Kbin = rand(1,L_Kbin)>0.5;
g = floor(L_Kbin/L_P_sum);
% 
% K= zeros(1,512);
for n = 0:L_P_sum-1
    K_n = [Kbin(g*n+1:g*n+g) P_sum_bin(n+1)];
    K1((g+1)*n+1:(g+1)*(n+1)) = K_n;
end 
K = [K1 Kbin(g*n+g+1:end)];
end

