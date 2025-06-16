function [qz_update,m_update,P_update] = CT_ekf_update(z,model,m,P)
plength= size(m,2);
zlength= size(z,2);

qz_update= zeros(plength,zlength);
m_update = zeros(model.x_dimCT,plength,zlength);
P_update = zeros(model.x_dimCT,model.x_dimCT,plength);

for idxp=1:plength
        [qz_temp,m_temp,P_temp] = ekf_update_single(z,model,m(:,idxp),P(:,:,idxp));
       qz_update(idxp,:)   = qz_temp;
        m_update(:,idxp,:) = m_temp;
        P_update(:,:,idxp) = P_temp;
end
end

function [qz_temp,m_temp,P_temp] = ekf_update_single(z,model,m,P)
%把omega删掉再更新xy，最后再加上
%可是协方差。。。
epsilon = 1e-10;
R_CT_omega =  [25  0    0;      %model.R_CT=[25 0; 0 25]多了一个omega数组不兼容
                0 25    0;
                0  0 pi/180];

H_noOmega = [1 0 0 0;   
             0 0 1 0];
m_noOmega = [m([1 2 3 4],:)];
P_noOmega = P(1:4, 1:4, :);

%% 5维求P
mu = model.H_CT*m;                  %预测观测5*n
S= R_CT_omega + model.H_CT*P*model.H_CT';                 %5*5*n
S= (S+ S')/2;   % addition step to avoid numerical problem
[S_positive] = makeMatrixPositiveDefinite(S, epsilon);
Vs= chol(S_positive); 
det_S= prod(diag(Vs))^2; 
inv_sqrt_S= inv(Vs); 
iS= inv_sqrt_S*inv_sqrt_S';
K  = P*model.H_CT'*iS;

P_temp = (eye(size(P))-K*model.H_CT)*P;%√


%% 4维求m
mu_noOmega = H_noOmega*m_noOmega;          %4*n
S_noOmega = model.R_CT + H_noOmega*P_noOmega*H_noOmega';  %4*4*n
S_noOmega= (S_noOmega+ S_noOmega')/2;
[S_positive_noOmega] = makeMatrixPositiveDefinite(S_noOmega, epsilon);
Vs_noOmega= chol(S_positive_noOmega); 
det_S_noOmega= prod(diag(Vs_noOmega))^2; 
inv_sqrt_S_noOmega= inv(Vs_noOmega); 
iS_noOmega= inv_sqrt_S_noOmega*inv_sqrt_S_noOmega';
K_noOmega  = P_noOmega*H_noOmega'*iS_noOmega;

qz_temp_noOmega = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S_noOmega) - 0.5*dot(z-repmat(mu_noOmega,[1 size(z,2)]),iS_noOmega*(z-repmat(mu_noOmega,[1 size(z,2)]))))';
%qz_temp好像就是似然啊

m_temp_noOmega = repmat(m_noOmega,[1 size(z,2)]) + K_noOmega*(z-repmat(mu_noOmega,[1 size(z,2)]));
% P_temp_withoutOmega = (eye(size(P))-K*model.H_CT)*P;

%% 
qz_temp = qz_temp_noOmega;%√
m_temp = [m_temp_noOmega; repmat(model.omega, 1, size(m_temp_noOmega, 2))];%√
% m_temp = [m_temp_withoutOmega;model.omega];
end



