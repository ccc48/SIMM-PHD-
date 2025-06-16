function z_gate= CT_gate_meas_ekf(z,gamma,model,m,P)
valid_idx = [];
zlength = size(z,2); if zlength==0, z_gate= []; return; end
plength = size(m,2); 

epsilon = 1e-10;
H_noOmega = [1 0 0 0;
             0 0 1 0];
P_noOmega = P(1:4, 1:4, :);
for j=1:plength
    %METHOD ONE
    Sj= model.R_CT + H_noOmega*P_noOmega(:,:,j)*H_noOmega';
    
    %METHOD TWO
    % [C_ekf,U_ekf]= CT_ekf_update_mat(model,m(:,j));
    % Sj= U_ekf*model.R_CT*U_ekf' + C_ekf*P(:,:,j)*C_ekf';

    Sj= (Sj+ Sj')/2;
    [Sj_positive] = makeMatrixPositiveDefinite(Sj, epsilon);
    Vs= chol(Sj_positive);
    det_Sj= prod(diag(Vs))^2; 
    inv_sqrt_Sj= inv(Vs);
    iSj= inv_sqrt_Sj*inv_sqrt_Sj'; 
    nu= z- repmat(CT_gen_observation_fn(model,m(:,j),zeros(size(model.D_CT,2),1)),[1 zlength]);
    dist= sum((inv_sqrt_Sj'*nu).^2);
    valid_idx= union(valid_idx,find( dist < gamma ));
end
z_gate = z(:,valid_idx);
end