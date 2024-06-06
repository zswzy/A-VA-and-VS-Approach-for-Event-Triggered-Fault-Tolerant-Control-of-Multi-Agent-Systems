function [Bfb, Cfb] = construct_fault_matrix(fault_config_act, fault_config_sens, B, C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


[N,m] = size(fault_config_act);
[~,r] = size(fault_config_sens);

Bfb = [];
Cfb = [];

for i = 1:N
    fault_config_act_i = fault_config_act(i,:);
    Bfi = B*diag(fault_config_act_i);
    Bfb = blkdiag(Bfb, Bfi);

    fault_config_sens_i = fault_config_sens(i,:);
    Cfi = diag(fault_config_sens_i)*C;
    Cfb = blkdiag(Cfb, Cfi);

end


end

