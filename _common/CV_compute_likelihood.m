function LK = CV_compute_likelihood(k,model, Z,est, m_predict, P_predict, w_predict)
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

    for l = 1:num_GM
        % 分量预测参数
        Hm = model.H_CV * m_predict(:, l);       % 预测量测均值
        S = model.H_CV * P_predict(:, :, l) * model.H_CV' + model.R_CV;  % 残差协方差
        S= (S+ S')/2;   % addition step to avoid numerical problem
        epsilon = 1e-10;
        [Sj_positive] = makeMatrixPositiveDefinite(S, epsilon);
        Vs= chol(Sj_positive);

        
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
    % LK1 = total_likelihood;


    % % 加入指数平滑（α=0.9，保留历史权重）
    % alpha = 0.9;
    % if k > 1
    %     LK = alpha * LK1 + (1-alpha) * est.Lk{k-1, 1};
    % else
    %     LK = LK1;
    % end
    
    % % 归一化模型概率
    % c_norm = sum(Lk_last .* c);
    % miu = (Lk_last .* c) / c_norm;
end

% function LK = CV_compute_likelihood(k,model,Z,est,m_predict,P_predict)

%     n = model.z_dim;
%     num_GM = size(m_predict,2);
% 
%     if isempty(m_predict)|isempty(Z)
%         LK= [];
%     else
% 
%     %残差  gating后，量测为空则残差为空
%     numZ = size(Z,2);
%     LK=[];
%     for Gll = 1:num_GM
%         v = Z -model.H*repmat(m_predict(:,Gll),[1 numZ]);
%         %残差为空则似然为空
%         %残差协方差
%         S = model.H *P_predict(:,:,Gll)*model.H' + model.R;
%         Vs= chol(S);        %Cholesky分解
%         det_S= prod(diag(Vs))^2; %S的行列式
%         inv_sqrt_S= inv(Vs);    %S的逆的平方根
%         iS= inv_sqrt_S*inv_sqrt_S';    %S的逆
%         %残差为空，则似然为空
%         for vll = 1:size(v,2)
%             Likel = exp(-0.5*n*log(2*pi)-0.5*log(det_S)-0.5*v(:,vll)'*iS*v(:,vll));  %似然
% 
%             %存
%             LK= [LK  Likel];
% 
%         end
%     end
%     end 
% end