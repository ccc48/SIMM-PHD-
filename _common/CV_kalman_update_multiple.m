function [qz_update,m_update,P_update] = CV_kalman_update_multiple(z,model,m,P)
plength= size(m,2);
zlength= size(z,2);

qz_update= zeros(plength,zlength);
m_update = zeros(model.x_dimCV,plength,zlength);
P_update = zeros(model.x_dimCV,model.x_dimCV,plength);

for idxp=1:plength
        [qz_temp,m_temp,P_temp] = kalman_update_single(z,model.H_CV,model.R_CV,m(:,idxp),P(:,:,idxp));
       qz_update(idxp,:)   = qz_temp;
        m_update(:,idxp,:) = m_temp;
        P_update(:,:,idxp) = P_temp;
end
end

function [qz_temp,m_temp,P_temp] = kalman_update_single(z,H,R,m,P)
epsilon = 1e-10;

mu = H*m;
S  = R+H*P*H'; 
S= (S+ S')/2;
[Sj_positive] = makeMatrixPositiveDefinite(S, epsilon);
Vs= chol(Sj_positive);
% Vs= chol(S); 
det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
K  = P*H'*iS;

qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*dot(z-repmat(mu,[1 size(z,2)]),iS*(z-repmat(mu,[1 size(z,2)]))))';
m_temp = repmat(m,[1 size(z,2)]) + K*(z-repmat(mu,[1 size(z,2)]));
P_temp = (eye(size(P))-K*H)*P;
end

