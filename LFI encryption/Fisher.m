function [P_new] = Fisher(S,N)
P_new = zeros(1,N);
P = linspace(1,N,N);
    for  flag = 1:length(S)
       P_new = P;
       P(flag) = P(S(flag));
       P(S(flag))= P_new(flag);    
    end
    P_new = P ;
end