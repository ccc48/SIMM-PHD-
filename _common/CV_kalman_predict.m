function [m_predict,P_predict] = CV_kalman_predict(model,m,P)      
plength= size(m,2);

m_predict = zeros(size(m));
P_predict = zeros(size(P));

for idxp=1:plength
    [m_temp,P_temp] = kalman_predict_single(model.F_CV,model.Q_CV,m(:,idxp),P(:,:,idxp));
    m_predict(:,idxp) = m_temp;
    P_predict(:,:,idxp) = P_temp;
end
end

function [m_predict,P_predict] = kalman_predict_single(F,Q,m,P)
m_predict = F*m;                % 状态预测
P_predict = Q + F*P*F';
end