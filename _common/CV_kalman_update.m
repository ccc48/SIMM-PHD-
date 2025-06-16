function [qz_update,m_update,P_update] = CV_kalman_update(z,model,m,P)

plength= size(m,2);    %预测状态：对应各高斯分量的个数
zlength= size(z,2);    %测量的个数

qz_update= zeros(plength,zlength);                 %初始化
m_update = zeros(model.x_dim,plength,zlength);
P_update = zeros(model.x_dim,model.x_dim,plength);

for idxp=1:plength
        [qz_temp,m_temp,P_temp] = kalman_update_single(z,model.H_CV,model.R_CV,m(:,idxp),P(:,:,idxp));
        qz_update(idxp,:)   = qz_temp;
        m_update(:,idxp,:)  = m_temp;
        P_update(:,:,idxp)  = P_temp;
end

function [qz_temp,m_temp,P_temp] = kalman_update_single(z,H,R,m,P)

 mu = H*m;           %预测观测
 S  = R+H*P*H';     %预测观测协方差
 Vs= chol(S);        %Cholesky分解
 det_S= prod(diag(Vs))^2; %S的行列式
 inv_sqrt_S= inv(Vs);    %S的逆的平方根
 iS= inv_sqrt_S*inv_sqrt_S';    %S的逆

K  = P*H'*iS;    %卡尔曼增益

% nu = z- repmat(mu,[1 size(z,2)])
% 
% nnnn = dot(z-repmat(mu,[1 size(z,2)]),iS*(z-repmat(mu,[1 size(z,2)])))
% AAA = iS*(z-repmat(mu,[1 size(z,2)]))
qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*dot(z-repmat(mu,[1 size(z,2)]),iS*(z-repmat(mu,[1 size(z,2)]))))';
%主要计算了观测值 z 在给定当前估计 m 时的概率密度值。
%1  高斯分布的常数部分：-0.5*size(z,1)*log(2*pi)    ----size(z,1)即z的行数/维度
%----%  计算多元高斯分布中与维度相关的常数部分。

%2  行列式的贡献：-0.5*log(det_S) 
%---%   通过 S 的行列式计算协方差矩阵 S 的复杂度或“广度”。这影响了概率密度的分布，行列式越大，密度越分散。

%3  dot(z-repmat(mu,[1 size(z,2)]),iS*(z-repmat(mu,[1 size(z,2)]))) 
%


m_temp = repmat(m,[1 size(z,2)]) + K*(z-repmat(mu,[1 size(z,2)]));


P_temp = (eye(size(P))-K*H)*P;     %更新协方差

