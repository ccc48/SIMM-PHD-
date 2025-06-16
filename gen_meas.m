function meas= gen_meas(model,truth)

%variables
meas.K= truth.K;
meas.Z= cell(truth.K,1);

% meas_noise = model.D_CV*randn(2,meas.K);%D_CV\CT\CA都一样[10 0;0 10]%%每时刻truth只有1个
% 计算距离(看全局噪声是否正态分布）
% x_noise = meas_noise(1, :);
% y_noise = meas_noise(2, :);
% distances = sqrt(x_noise.^2 + y_noise.^2);

% meas_noise = zeros(2,meas.K);%零噪声

%generate measurements
for k=1:truth.K
    if truth.N(k) > 0
        idx= find( rand(truth.N(k),1) <= model.P_D );%生成0-1之间随机数，判断目标是否被检测到                                            %detected target indices
        meas.Z{k}= gen_observation_fn(model,truth.X{k}(:,idx),k);                          %single target observations if detected 
    end
    N_c= poissrnd(model.lambda_c);                                                               %number of clutter points
    C= repmat(model.range_c(:,1),[1 N_c])+ diag(model.range_c*[ -1; 1 ])*rand(model.z_dim,N_c);  %clutter generation
    meas.Z{k}= [ meas.Z{k} C ];                                                                  %measurement is union of detections and clutter
end
end  