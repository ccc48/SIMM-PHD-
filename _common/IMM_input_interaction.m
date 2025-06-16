function [w01,m01,P01,c]=IMM_input_interaction(k,model,est,w,m,P,L,miu)
    mu = zeros(model.M,model.M);    %混合概率矩阵
    c = zeros(model.M,1);     %归一化常数C_j

    for j = 1:model.M
        c(j,:) =  sum( model.Pmn(:,j).*miu );
    end

    for i = 1:model.M
        for j = 1:model.M
            mu(i,j) = model.Pmn(i,j) *miu(i)./c(j);
        end
    end

    %
    if k==1
        w_1M = w;  w_2M = w;  
        m_1M = m;  m_2M = m;  
        P_1M = P;  P_2M = P; 
        N_1M = L;  N_2M = L; 
        NN = [N_1M  N_2M ];
    else
        w_1M = est.W{k-1,1};  w_2M = est.W{k-1,2}; 
        m_1M = est.X{k-1,1};  m_2M = est.X{k-1,2};  
        P_1M = est.P{k-1,1};  P_2M = est.P{k-1,2};  
        N_1M = est.L{k-1,1};  N_2M = est.L{k-1,2}; 
    end
  
    % 计算每个CV分量与每个CT分量和每个CA分量的交互输入
    for j = 1:model.M
        w01{j} = w_1M.*mu(1,j) + w_2M.*mu(2,j) ;
        m01{j} = m_1M.*mu(1,j) + m_2M.*mu(2,j) ;
    end
    for j = 1:model.M
        P01{j} = mu(1,j) * (P_1M + (m_1M - m01{j}) * (m_1M - m01{j})')... 
                +mu(2,j) * (P_2M + (m_2M - m01{j}) * (m_2M - m01{j})');
            
    end
        
end