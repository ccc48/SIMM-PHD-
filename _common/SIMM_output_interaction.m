function [w_update,m_update,P_update,L_update] = SIMM_output_interaction(k,model,estS,miu)

w = estS.W;  %M个模型的wmPL
m = estS.X;
P = estS.P;
L = estS.L;

w_update = [];    %交互输出
m_update = [];
P_update = [];
L_update = [];

N = cell2mat(L(k,1));%高斯分量个数

% 计算每个CV分量与每个CT分量和每个CA分量的交互输出
[~,jmax] = max(miu);  %直接找出可能性最大的model
w_update = w{k,jmax};      %可能性最大的model就是对应的交互输出
m_update = m{k,jmax};      %wmP
P_update = P{k,jmax};
   
L_update = length(w_update);
end