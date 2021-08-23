function [P_old] = RFisher(S,N)
S_n = zeros(1,N);
P_old = zeros(1,N);
P = linspace(1,N,N);
for  flag = 1:length(S)
       P_new = P;
       P(flag) = P(S(flag));
       P(S(flag))= P_new(flag);    
    end
    P_new = P ;
    for i = 1:length(S)
        for S_flag = 1:length(S)
            if P_new(i) == S_flag
                S_n(S_flag)=i;
                S_flag=S_flag+1;
            else
                 S_flag=S_flag+1;
            end
        end           
    end
    P_old = S_n;