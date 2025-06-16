function z_gate= CV_gate_meas_gms(z,gamma,model,m,P)

valid_idx = [];
zlength = size(z,2); if zlength==0, z_gate= []; return; end
plength = size(m,2);

for j=1:plength
    Sj= model.R_CV + model.H_CV*P(:,:,j)*model.H_CV';
    Sj= (Sj+ Sj')/2;   % addition step to avoid numerical problem
 
    epsilon = 1e-10;
    [S_positive] = makeMatrixPositiveDefinite(Sj, epsilon);
    Vs= chol(S_positive); 
    % Vs= chol(Sj); 
    det_Sj= prod(diag(Vs))^2; inv_sqrt_Sj= inv(Vs);
    iSj= inv_sqrt_Sj*inv_sqrt_Sj';
    nu= z- model.H_CV*repmat(m(:,j),[1 zlength]);
    dist= sum((inv_sqrt_Sj'*nu).^2);
    valid_idx= union(valid_idx,find( dist < gamma ));
end
z_gate = z(:,valid_idx);
end