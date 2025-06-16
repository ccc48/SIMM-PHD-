function [w01,m01,P01,c]=SIMM_input_interaction(k,model,estS,w,m,P,L,miu)
    mu = zeros(model.M,model.M);    % 混合概率矩阵
    c = zeros(model.M,1);     % 归一化常数C_j

    % 计算归一化常数 (使用max而非IMM中的sum)
    for j = 1:model.M
        c(j,:) = max(model.Sij(:,j).*miu);
    end

    % 计算混合概率
    for i = 1:model.M
        for j = 1:model.M
            mu(i,j) = model.Sij(i,j) * miu(i) / c(j);
        end
    end

    % 获取前一时刻各模型的状态
    if k==1
        % 初始时刻使用相同的初始值
        w_1M = w;  w_2M = w;  
        m_1M = m;  m_2M = m;  
        P_1M = P;  P_2M = P;  
    else
        % 非初始时刻使用前一时刻估计值
        w_1M = estS.W{k-1,1};  w_2M = estS.W{k-1,2};  
        m_1M = estS.X{k-1,1};  m_2M = estS.X{k-1,2}; 
        P_1M = estS.P{k-1,1};  P_2M = estS.P{k-1,2}; 
    end
    WW = {w_1M  w_2M };
    MM = {m_1M  m_2M };
    PP = {P_1M  P_2M };

    % 计算交互输入 - 加权混合所有模型的状态
    for j = 1:model.M
        [~,lmax] = max(mu(:,j));  %
        w01{j} = WW{lmax};        %权重、状态
        m01{j} = MM{lmax};
    end
    for j = 1:model.M    
        [~,l] = max(mu(:,j));  
        P01{j} = PP{j} + [MM{j} - m01{l}]*[MM{j} - m01{l}]';%协方差
    end
        
end