function [I_RGB] = decryption(K,I_R,I_G,I_B,XR_xy,XG_xy,XB_xy,V1,V2,V3)
[M,N]=size(I_R);
L = M*N;
[x1,p1,x2,p2] = generate_x0_p(K);
 x0 = mod(x1+x2,1);
 p = mod(p1+p2,0.5);
 %% 生成加密所需的伪随机序列
L_c=N/2;
L_l=M/2; 
L_S=L_c+L_l+1000;
S_c = zeros(1,L_S);
S_c(1) = pwlcm(x1,p1);
 for n=2:L_S
     S_c(n) = pwlcm(S_c(n-1),p1);
 end
S1_c1 = mod(floor(S_c*power(10,15)),L_l)+1;
S1_c2 = mod(floor(S_c*power(10,15)),L_c)+1;
S_J_M = S1_c1(1001:L_l+1000);
S_J_N = S1_c2(L_l+1001:L_l+L_c+1000);
L_D=200+L_c*3 + L_l; 
S_d = zeros(1,L_D); 
S_d(1) = pwlcm(x2,p2);
 for n=2:L_D
     S_d(n) = pwlcm(S_d(n-1),p2);
 end
A1 = mod(floor(S_d*power(10,15)),L_l)+1; 
A2 = mod(floor(S_d*power(10,15)),L_c*3)+1; 
A = A1(201:end); 
A_l =A(1:L_l); 
A_c = A2(L_l+201:end);
L_K=M*N+L_c*L_l+200;
S_k = zeros(1,L_K);
S_k(1) = pwlcm(x0,p);
 for n=2:L_K
     S_k(n) = pwlcm(S_k(n-1),p);
 end 
S_k_LL = S_k(201:L_c*L_l+200);
SS = mod(floor(S_k_LL*power(10,15)),510); 
S_others = S_k(L_c*L_l+201:end);
SS_s = mod(floor(S_others*power(10,15)),256);
%% 解密
I_R=double(I_R);
I_G=double(I_G);
I_B=double(I_B);
DDR_s = reshape(I_R,1,[]);
ER_s = mod(255*2+DDR_s-V1-SS_s,255);
for n = M*N-1:-1:1
    ER_s(n) = mod(255*2+DDR_s(n)-SS_s(n)-DDR_s(n+1),255);
end
FR_s = zeros(1,M*N);
FR_s(1) = mod(255*2+ER_s(1)-V1-SS_s(1),255);
for n = 2:M*N
    FR_s(n) = mod(255*2+ER_s(n)-SS_s(n)-ER_s(n-1),255); 
end
DDG_s = reshape(I_G,1,[]);
EG_s = mod(255*2+DDG_s-V2-SS_s,255);
for n = M*N-1:-1:1
    EG_s(n) = mod(255*2+DDG_s(n)-SS_s(n)-DDG_s(n+1),255);
end
FG_s = zeros(1,M*N);
FG_s(1) = mod(255*2+EG_s(1)-V2-SS_s(1),255);
for n = 2:M*N
    FG_s(n) = mod(255*2+EG_s(n)-SS_s(n)-EG_s(n-1),255); 
end
DDB_s = reshape(I_B,1,[]);
EB_s = mod(255*2+DDB_s-V3-SS_s,255);
for n = M*N-1:-1:1
    EB_s(n) = mod(255*2+DDB_s(n)-SS_s(n)-DDB_s(n+1),255);
end
FB_s = zeros(1,M*N);
FB_s(1) = mod(255*2+EB_s(1)-V3-SS_s(1),255);
for n = 2:M*N
    FB_s(n) = mod(255*2+EB_s(n)-SS_s(n)-EB_s(n-1),255); 
end
X_R = zeros(1,L);
for n = 1:L
    if XR_xy(n) == 0
        X_R(n)  = mod(FR_s(n),255);
    else
        if XR_xy(n) < 0
           X_R(n) = FR_s(n) + 255;
        else
          X_R(n) = FR_s(n) - 255;   
        end      
    end
end
FFR_s = reshape(X_R,[M,N]);
X_G = zeros(1,L);
for n = 1:L
    if XG_xy(n) == 0
        X_G(n)  = mod(FG_s(n),255);
    else
        if XG_xy(n) < 0
           X_G(n) = FG_s(n) + 255;
        else
          X_G(n) = FG_s(n) - 255;   
        end      
    end
end
FFG_s = reshape(X_G,[M,N]);
X_B = zeros(1,L);
for n = 1:L
    if XB_xy(n) == 0
        X_B(n)  = mod(FB_s(n),255);
    else
        if XB_xy(n) < 0
           X_B(n) = FB_s(n) + 255;
        else
          X_B(n) = FB_s(n) - 255;   
        end      
    end
end
FFB_s = reshape(X_B,[M,N]);
[LLR4,LHR4,HLR4,HHR4]= dwt2(FFR_s,'haar');
[LLG4,LHG4,HLG4,HHG4]= dwt2(FFG_s,'haar');
[LLB4,LHB4,HLB4,HHB4]= dwt2(FFB_s,'haar');
imgR = [LHR4 HLR4 HHR4];
imgG = [LHG4 HLG4 HHG4];
imgB = [LHB4 HLB4 HHB4];
X_l_c = RFisher(A_l,L_l);
X_c_c = RFisher(A_c,L_c*3);
I_NR = scramble_N(imgR,X_c_c);
I_MR = scramble_M(I_NR,X_l_c);
LHR3 = I_MR(:,1:N/2);
HLR3 = I_MR(:,N/2+1:N);
HHR3 = I_MR(:,N+1:end);
I_NG = scramble_N(imgG,X_c_c);
I_MG = scramble_M(I_NG,X_l_c);
LHG3 = I_MG(:,1:N/2);
HLG3 = I_MG(:,N/2+1:N);
HHG3 = I_MG(:,N+1:end);
I_NB = scramble_N(imgB,X_c_c);
I_MB= scramble_M(I_NB,X_l_c);
LHB3 = I_MB(:,1:N/2);
HLB3 = I_MB(:,N/2+1:N);
HHB3 = I_MB(:,N+1:end);
X_N_c = RFisher(S_J_N,L_c);
X_M_c = RFisher(S_J_M,L_l);
I_RN = scramble_N(LLR4,X_N_c);
I_LLR3 = scramble_M(I_RN,X_M_c);
I_GN = scramble_N(LLG4,X_N_c);
I_LLG3 = scramble_M(I_GN,X_M_c);
I_BN = scramble_N(LLB4,X_N_c);
I_LLB3 = scramble_M(I_BN,X_M_c);
LLR3 = reshape(I_LLR3,1,[]);
LLG3 = reshape(I_LLG3,1,[]);
LLB3 = reshape(I_LLB3,1,[]);
ER = mod(511*2+LLR3-V1-SS,511);
for n = L_l*L_c-1:-1:1
    ER(n) = mod(511*2+LLR3(n)-SS(n)-LLR3(n+1),511);
end
FR = zeros(1,L_l*L_c);
FR(1) = mod(511*2+ER(1)-V1-SS(1),511);
for n = 2:L_l*L_c
    FR(n) = mod(511*2+ER(n)-SS(n)-ER(n-1),511); 
end
FFR = reshape(FR,[L_l,L_c]);
EG = mod(510*2+LLG3-V2-SS,510);
for n = L_l*L_c-1:-1:1
    EG(n) = mod(510*2+LLG3(n)-SS(n)-LLG3(n+1),510);
end
FG = zeros(1,L_l*L_c);
FG(1) = mod(510*2+EG(1)-V2-SS(1),510);
for n = 2:L_l*L_c
    FG(n) = mod(510*2+EG(n)-SS(n)-EG(n-1),510); 
end
FFG = reshape(FG,[L_l,L_c]);
EB = mod(510*2+LLB3-V3-SS,510);
for n = L_l*L_c-1:-1:1
    EB(n) = mod(510*2+LLB3(n)-SS(n)-LLB3(n+1),510);
end
FB = zeros(1,L_l*L_c);
FB(1) = mod(510*2+EB(1)-V3-SS(1),510);
for n = 2:L_l*L_c
    FB(n) = mod(510*2+EB(n)-SS(n)-EB(n-1),510); 
end
FFB = reshape(FB,[L_l,L_c]);
X_N_c = RFisher(S_J_N,L_c);
X_M_c = RFisher(S_J_M,L_l);
IR_N = scramble_N(FFR,X_N_c);
LLR_3 = scramble_M(IR_N,X_M_c);
IG_N = scramble_N(FFG,X_N_c);
LLG_3 = scramble_M(IG_N,X_M_c);
IB_N = scramble_N(FFB,X_N_c);
LLB_3 = scramble_M(IB_N,X_M_c);
IR= idwt2(LLR_3,LHR3,HLR3,HHR3,'haar');
IG= idwt2(LLG_3,LHG3,HLG3,HHG3,'haar');
IB= idwt2(LLB_3,LHB3,HLB3,HHB3,'haar');
I_RGB(:,:,1)=uint8(IR);
I_RGB(:,:,2)=uint8(IG);
I_RGB(:,:,3)=uint8(IB);
end
