function [w_update,m_update,P_update,L_update] = IMM_output_interaction(k,model,est,miu)

w = est.W;
m = est.X;
P = est.P;
L = est.L;

w_update = [];
m_update = [];
P_update = [];
L_update = [];

N = cell2mat(L(k,1));%高斯分量个数

% 计算每个CV分量与每个CT分量和每个CA分量的交互输出
for mll = 1:N 
    %交互
    w_j = w{k,1}(mll)*miu(1) + w{k,2}(mll)*miu(2) + w{k,3}(mll)*miu(3);
    m_j = m{k,1}(:,mll)*miu(1) + m{k,2}(:,mll)*miu(2) + m{k,3}(:,mll)*miu(3);
    P_j = miu(1)*(P{k,1}(:,:,mll) +(m{k,1}(:,mll) - m_j)*(m{k,1}(:,mll) - m_j)')...
        + miu(2)*(P{k,2}(:,:,mll) +(m{k,2}(:,mll) - m_j)*(m{k,2}(:,mll) - m_j)')...
        + miu(3)*(P{k,3}(:,:,mll) +(m{k,3}(:,mll) - m_j)*(m{k,3}(:,mll) - m_j)') ;
    %存
    w_update = [w_update; w_j]; 
    m_update = [m_update m_j];
    P_update = cat(3,P_update, P_j);
end
%看最后有几个
L_update = length(w_update);
end


