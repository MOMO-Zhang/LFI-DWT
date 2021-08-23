function [x1,p1,x2,p2] = generate_x0_p(K)
K_d = zeros(1,128);
for n = 0:127
    K_d(n+1) = bin2dec(int2str(K(4*n+1:4*n+4))); 
end


K_d_mean = sum(K_d)/power(2,11);
 K_32=0;
 for n = 1:32
     K_32 = K_32+K_d(n);
 end
 
x1 = mod((K_d_mean + K_32/480 + 0.001),1);

K_64=0;
 for n = 34:64
     K_64 = K_64 + K_d(n);
 end
p1= mod(K_d_mean + K_64/480 + x1,0.5);

K_96=0;
 for n = 66:96
     K_96 = K_96 + K_d(n);
 end
x2 = mod((K_d_mean + K_96/480 + p1),1);

K_128=0;
 for n = 98:128
     K_128 = K_128 + K_d(n);
 end
p2 = mod(K_d_mean + K_128/480 + x2,0.5);
 for n = 1:100
       x1 = pwlcm(x1,p1);
 end
 
 for n = 1:100
       x2 = pwlcm(x2,p2);
 end
 
%  x0 = mod(x1+x2,1);
%  p = mod(p1+p2,0.5);
end

