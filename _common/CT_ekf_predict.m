function [m_predict,P_predict] = CT_ekf_predict(model,m,P)    
plength= size(m,2);
%CT，目标状态4——>5
% 创建一个新行，每个元素都是omega
new_row = model.omega * ones(1, plength);  % 1×n的行向量
% 合并提取的行和新行
X_1i = [m; new_row];  % 结果为5×n的矩阵

P_1i = zeros(5, 5, plength);% 创建全零的5×5×n矩阵
P_1i(1:4, 1:4, :) = P;% 复制原始矩阵到新矩阵的左上角

% 设置每个切片的新增对角线元素（第5行第5列）
for i = 1:plength
    P_1i(5, 5, i) = model.omega;  % 用你需要的值替换desired_value
end
m_predict = zeros(size(X_1i));
P_predict = zeros(size(P_1i));

for idxp=1:plength
    [m_temp,P_temp] = ekf_predict_single(model,X_1i(:,idxp),P_1i(:,:,idxp));
    m_predict(:,idxp) = m_temp;
    P_predict(:,:,idxp) = P_temp;
end
end

function [m_predict,P_predict] = ekf_predict_single(model,m,P)
%METHOD ONE 状态转移矩阵
m_predict = model.F_CT * m;
P_predict= model.G_CT*model.Q_CT*model.G_CT' + model.F_CT*P*model.F_CT';
P_predict= (P_predict+P_predict')/2;                    % addition step to avoid numerical problem

%METHOD TWO 变化的状态转移矩阵
% m_predict = CT_gen_newstate_fn(model,m,'noiseless');
% [F_ekf,G_ekf]= CT_ekf_predict_mat(model,m);                % user specified function for application
% P_predict= G_ekf*model.Q_CT*G_ekf' + F_ekf*P*F_ekf';
% P_predict= (P_predict+P_predict')/2;                    % addition step to avoid numerical problem
end