function [x] = pwlcm(x,p)
% PWLCM发生器
%   此处显示详细说明
if x > 0&&x <= p
    x= x/p;
elseif x > p&&x <= 0.5
    x = (x-p)/(0.5-p);
else
    x = pwlcm(1-x,p);
end