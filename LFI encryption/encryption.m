function [I_R,I_G,I_B,XR_xy,XG_xy,XB_xy,V1,V2,V3] = encryption(K,I)
I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B
[M,N]=size(I1);
L = M*N;
P1_sum=sum(I1(:));
V1=mod(P1_sum,256);
P2_sum=sum(I2(:));
V2=mod(P2_sum,256);
P3_sum=sum(I3(:));
V3=mod(P3_sum,256);
[x1,p1,x2,p2] = generate_x0_p(K);
 x0 = mod(x1+x2,1);
 p = mod(p1+p2,0.5);
%% S
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
%%dwt
[LLR,LHR,HLR,HHR]= dwt2(I1,'haar');
[LLG,LHG,HLG,HHG]= dwt2(I2,'haar');
[LLB,LHB,HLB,HHB]= dwt2(I3,'haar');
%%
X_M = Fisher(S_J_M,L_l);
X_N = Fisher(S_J_N,L_c);
LLR_M = scramble_M(LLR,X_M);
LLR_fisher = scramble_N(LLR_M,X_N);
LLR1 = reshape(LLR_fisher,[],1);
LLG_M = scramble_M(LLG,X_M);
LLG_fisher = scramble_N(LLG_M,X_N);
LLG1 = reshape(LLG_fisher,[],1);
LLB_M = scramble_M(LLB,X_M);
LLB_fisher = scramble_N(LLB_M,X_N);
LLB1 = reshape(LLB_fisher,[],1);
C_R = zeros(1,L_l*L_c);
C_R(1) = mod(V1+SS(1)+LLR1(1),511);
for n = 2:L_l*L_c
    C_R(n) = mod(LLR1(n)+SS(n)+C_R(n-1),511); 
end
D_R = mod(V1+SS+C_R,511);
for n = L_l*L_c-1:-1:1
    D_R(n) = mod(D_R(n+1)+SS(n)+C_R(n),511);
end
DD_R = reshape(D_R,[L_l,L_c]);

C_G = zeros(1,L_l*L_c);
C_G(1) = mod(V2+SS(1)+LLG1(1),510);
for n = 2:L_l*L_c
    C_G(n) = mod(LLG1(n)+SS(n)+C_G(n-1),510); 
end
D_G = mod(V2+SS+C_G,510);
for n = L_l*L_c-1:-1:1
    D_G(n) = mod(D_G(n+1)+SS(n)+C_G(n),510);
end
DD_G = reshape(D_G,[L_l,L_c]);

C_B = zeros(1,L_l*L_c);
C_B(1) = mod(V3+SS(1)+LLB1(1),510);
for n = 2:L_l*L_c
    C_B(n) = mod(LLB1(n)+SS(n)+C_B(n-1),510); 
end
D_B = mod(V3+SS+C_B,510);
for n = L_l*L_c-1:-1:1
    D_B(n) = mod(D_B(n+1)+SS(n)+C_B(n),510);
end
DD_B = reshape(D_B,[L_l,L_c]);
X1_M = Fisher(S_J_M,L_l);
X1_N = Fisher(S_J_N,L_c);
LLR1_M = scramble_M(DD_R,X1_M);
LLR_1 = scramble_N(LLR1_M,X1_N);
LLG1_M = scramble_M(DD_G,X1_M);
LLG_1 = scramble_N(LLG1_M,X1_N);
LLB1_M = scramble_M(DD_B,X1_M);
LLB_1 = scramble_N(LLB1_M,X1_N);
LLR_others = [LHR HLR HHR]; 
LLG_others = [LHG HLG HHG]; 
LLB_others = [LHB HLB HHB];
X_l_M = Fisher(A_l,L_l);
X_c_N = Fisher(A_c,L_c*3);
LLR_others_M = scramble_M(LLR_others,X_l_M);
LLR_others_N = scramble_N(LLR_others_M ,X_c_N);
LHR1 = LLR_others_N(:,1:N/2);
HLR1 = LLR_others_N(:,N/2+1:N);
HHR1 = LLR_others_N(:,N+1:end);
JMR_P= idwt2(LLR_1,LHR1,HLR1,HHR1,'haar');

LLG_others_M = scramble_M(LLG_others,X_l_M);
LLG_others_N = scramble_N(LLG_others_M ,X_c_N);
LHG1 = LLG_others_N(:,1:N/2);
HLG1 = LLG_others_N(:,N/2+1:N);
HHG1 = LLG_others_N(:,N+1:end);
JMG_P= idwt2(LLG_1,LHG1,HLG1,HHG1,'haar');

LLB_others_M = scramble_M(LLB_others,X_l_M);
LLB_others_N = scramble_N(LLB_others_M ,X_c_N)
LHB1 = LLB_others_N(:,1:N/2);
HLB1 = LLB_others_N(:,N/2+1:N);
HHB1 = LLB_others_N(:,N+1:end);
JMB_P= idwt2(LLB_1,LHB1,HLB1,HHB1,'haar');

JMR_S = reshape(JMR_P,[],1);
JMG_S = reshape(JMG_P,[],1);
JMB_S = reshape(JMB_P,[],1);

XR_xy = zeros(1,L);
for n = 1:L
    if JMR_S(n) < 0
        XR_xy(n) = 1;
    else
        if JMR_S(n)>255
            XR_xy(n) = -1;
        else
            XR_xy(n) = 0;
        end
    end
end
XG_xy = zeros(1,L);
for n = 1:L
    if JMG_S(n) < 0
        XG_xy(n) = 1;
    else
        if JMG_S(n)>255
            XG_xy(n) = -1;
        else
            XG_xy(n) = 0;
        end
    end
end
XB_xy = zeros(1,L);
for n = 1:L
    if JMB_S(n) < 0
        XB_xy(n) = 1;
    else
        if JMB_S(n)>255
            XB_xy(n) = -1;
        else
            XB_xy(n) = 0;
        end
    end
end
CR_s = zeros(1,M*N);
CR_s(1) = mod(V1+SS_s(1)+JMR_S(1),255);
for n = 2:M*N
    CR_s(n) = mod(JMR_S(n)+SS_s(n)+CR_s(n-1),255); 
end
DR_s = mod(V1+SS_s+CR_s,255);
for n = M*N-1:-1:1
    DR_s(n) = mod(DR_s(n+1)+SS_s(n)+CR_s(n),255);
end
I_R = reshape(DR_s,[M,N]);
CB_s = zeros(1,M*N);
CB_s(1) = mod(V3+SS_s(1)+JMB_S(1),255);
for n = 2:M*N
    CB_s(n) = mod(JMB_S(n)+SS_s(n)+CB_s(n-1),255); 
end
DB_s = mod(V3+SS_s+CB_s,255);
for n = M*N-1:-1:1
    DB_s(n) = mod(DB_s(n+1)+SS_s(n)+CB_s(n),255);
end
I_B = reshape(DB_s,[M,N]);
CG_s = zeros(1,M*N);
CG_s(1) = mod(V2+SS_s(1)+JMG_S(1),255);
for n = 2:M*N
    CG_s(n) = mod(JMG_S(n)+SS_s(n)+CG_s(n-1),255); 
end
DG_s = mod(V2+SS_s+CG_s,255);
for n = M*N-1:-1:1
    DG_s(n) = mod(DG_s(n+1)+SS_s(n)+CG_s(n),255);
end
I_G = reshape(DG_s,[M,N]);
end