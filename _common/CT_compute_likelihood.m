function LK = CT_compute_likelihood(k,model, Z,est, m_predict, P_predict, w_predict)
    % 输入参数：
    % model: 包含H, R, z_dim等模型参数
    % Z: 当前时刻门控后的量测集合 [z_dim x numZ]
    % m_predict: 预测的高斯分量均值 [x_dim x num_GM]
    % P_predict: 预测的高斯分量协方差 [x_dim x x_dim x num_GM]
    % omega_predict: 高斯分量的权重 [1 x num_GM]

    if isempty(m_predict) || isempty(Z)%其中任何一个为空，就把 LK 设为 0 并结束函数
        LK = 0;
        return;
    end

    n = model.z_dim;                 %量测的维度
    num_GM = size(m_predict, 2);     %高斯分量的数量
    numZ = size(Z, 2);     %当前时刻门控后的量测数量
    total_likelihood = 0;  %初始化为 0，用于后续累加每个高斯分量的似然和

    H_noOmega = [1 0 0 0;   
                 0 0 1 0];
    m_noOmega = [m_predict([1 2 3 4],:)];
    P_noOmega = P_predict(1:4, 1:4, :);

    for l = 1:num_GM
        Hm = H_noOmega * m_noOmega(:, l);       % 预测量测均值
        S = H_noOmega * P_noOmega(:, :, l) * H_noOmega' + model.R_CT;  % 残差协方差

        % eta = CT_gen_observation_fn(model,m_predict(:, l),'noiseless');
        % [H_ekf,U_ekf]= CT_ekf_update_mat(model,m_predict(:, l));                 % user specified function for application
        % S= U_ekf*model.R_CT*U_ekf'+H_ekf*P_predict(:, :, l)*H_ekf'; 
        S= (S+ S')/2;   % addition step to avoid numerical problem
        epsilon = 1e-10;
        [S_positive] = makeMatrixPositiveDefinite(S, epsilon);
        Vs= chol(S_positive); 
        % Vs = chol(S);                         % Cholesky分解,得到上三角矩阵 Vs
        det_S = prod(diag(Vs))^2;             % S 的行列式
        inv_S = inv(S);                       % 直接计算S 的逆矩阵
        
        log_likelihoods = zeros(1, numZ);      %第 l 个高斯分量对每个量测的似然
        for z = 1:numZ
            v = Z(:, z) - Hm;                  % 残差
            quad_form = v' * inv_S * v;        % 残差 v 的马氏距离
            log_likelihood = -0.5*(n*log(2*pi) + log(det_S) + quad_form);  %第z个量测在第l个高斯分量下的对数似然
            log_likelihoods(z) = exp(log_likelihood);%对对数似然进行指数运算，得到第 z 个量测在第 l 个高斯分量下的似然
        end
        
        % 加权求和：高斯分量权重 * 量测似然和
        total_likelihood = total_likelihood + w_predict(l) * sum(log_likelihoods);
    end

    LK = total_likelihood;
end