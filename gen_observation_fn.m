function Z = gen_observation_fn(model,x,k) 

    meas_noise = model.D_CV*randn(2,1);%D_CV\CT\CA都一样[10 0;0 10]%%每时刻truth只有1个
    Z1 = x([1 3],:);
    % Z = Z1 + meas_noise(:,k);
    Z = Z1 + meas_noise;

end