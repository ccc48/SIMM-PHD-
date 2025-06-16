function [w_update01,m_update01,  P_update01,  L_update01,LK01,...
          w_update02,m_update002, P_update002, L_update02,LK02,...
          w_update03,m_update003, P_update003, L_update03,LK03] = CVCTfilters(k,model,meas,est,filter,w0,m0,P0,L0)
          %CT  5-->4    
%=== Filtering 
%取出
%CV
w01 = w0{1};m01 = m0{1};P01 = P0{1};L01 = L0{1};
%CT
w02 = w0{2};m02 = m0{2};P02 = P0{2};L02 = L0{2};
%CT2
w03 = w0{3};m03 = m0{3};P03 = P0{3};L03 = L0{3};

%% ---prediction
%CV
[m_predict01,P_predict01] = CV_kalman_predict(model,m01,P01);       %surviving components
w_predict01= model.P_S*w01;                                                                  %surviving weights
w_predict01= cat(1,model.w_birth,w_predict01);                                                   %append birth weights
m_predict01= cat(2,model.m_birthCV,m_predict01); 
P_predict01= cat(3,model.P_birthCV,P_predict01);            %append birth components                                                 
L_predict01= model.L_birth+L01;                                                              %number of predicted components

%CT 
[m_predict02,P_predict02] = CT_ekf_predict(model,m02,P02);                          %surviving components
w_predict02= model.P_S*w02;                                                                  %surviving weights
w_predict02= cat(1,model.w_birth,w_predict02);
m_predict02= cat(2,model.m_birthCT,m_predict02); 
P_predict02= cat(3,model.P_birthCT,P_predict02);            %append birth components                                                      %append birth weights
L_predict02= model.L_birth+L02;

%CT 2
[m_predict03,P_predict03] = CT_ekf_predict(model,m03,P03);                          %surviving components
w_predict03= model.P_S*w03;                                                                  %surviving weights
w_predict03= cat(1,model.w_birth,w_predict03);
m_predict03= cat(2,model.m_birthCT2,m_predict03); 
P_predict03= cat(3,model.P_birthCT2,P_predict03);            %append birth components                                                      %append birth weights
L_predict03= model.L_birth+L03;

%% ---gating
%CV
if filter.gate_flag
    Z01 = CV_gate_meas_gms(meas.Z{k},filter.gamma,model,m_predict01,P_predict01);        
end
 
%CT
if filter.gate_flag
    Z02 = CT_gate_meas_ekf(meas.Z{k},filter.gamma,model,m_predict02,P_predict02);        
end

%CT 2
if filter.gate_flag
    Z03 = CT_gate_meas_ekf(meas.Z{k},filter.gamma,model,m_predict03,P_predict03);        
end

%% ---update
%number of measurements
m1= size(Z01,2);
m2= size(Z02,2);
m3= size(Z03,2);

%missed detection term漏检
%CV
w_update01 = model.Q_D*w_predict01;
m_update01 = m_predict01;
P_update01 = P_predict01;
L_update01 = L_predict01;

%CT
w_update02 = model.Q_D*w_predict02;
m_update02 = m_predict02;
P_update02 = P_predict02;
L_update02 = L_predict02;

%CT 2
w_update03 = model.Q_D*w_predict03;
m_update03 = m_predict03;
P_update03 = P_predict03;
L_update03 = L_predict03;

%检测
if m1~=0
%m detection terms 
%CV
LK01 = 0;
[qz_temp01,m_temp01,P_temp01] = CV_kalman_update_multiple(Z01,model,m_predict01,P_predict01);
    for ell=1:m1
        w_temp01 = model.P_D*w_predict01(:).*qz_temp01(:,ell); 
        w_temp01 = w_temp01./(model.lambda_c*model.pdf_c + sum(w_temp01));
        w_update01 = cat(1,w_update01,w_temp01);
        m_update01 = cat(2,m_update01,m_temp01(:,:,ell));
        P_update01 = cat(3,P_update01,P_temp01);
        L_update01 = length(w_update01);
        %%Likelihood
        %%%CV
        LK01 = LK01 + sum(qz_temp01(:,ell).*w_temp01);
    end
else
    %%Likelihood
    %%%CV
    LK01 = eps;
end

if m2~=0
%CT
LK02=0;
[qz_temp02,m_temp02,P_temp02] = CT_ekf_update(Z02,model,m_predict02,P_predict02);
    for ell=1:m2
        w_temp02 = model.P_D*w_predict02(:).*qz_temp02(:,ell); 
        w_temp02 = w_temp02./(model.lambda_c*model.pdf_c + sum(w_temp02));
        w_update02 = cat(1,w_update02,w_temp02);
        m_update02 = cat(2,m_update02,m_temp02(:,:,ell));
        P_update02 = cat(3,P_update02,P_temp02);
        L_update02 = length(w_update02);
        %%Likelihood
        %%%CT
        LK02 = LK02 + sum(qz_temp02(:,ell).*w_temp02);
    end
else
    %%Likelihood
    %%%CT
    LK02 = eps;
end         
            

if m3~=0
%CT 2
LK03=0;
[qz_temp03,m_temp03,P_temp03] = CT_ekf_update(Z03,model,m_predict03,P_predict03);
    for ell=1:m3
        w_temp03 = model.P_D*w_predict03(:).*qz_temp03(:,ell); 
        w_temp03 = w_temp03./(model.lambda_c*model.pdf_c + sum(w_temp03));
        w_update03 = cat(1,w_update03,w_temp03);
        m_update03 = cat(2,m_update03,m_temp03(:,:,ell));
        P_update03 = cat(3,P_update03,P_temp03);
        L_update03 = length(w_update03);
        %%Likelihood
        %%%CT
        LK03 = LK03 + sum(qz_temp03(:,ell).*w_temp03);
    end
else
    %%Likelihood
    %%%CT
    LK03 = eps;
end      
 
%% ---mixture management
%CV
L_posterior01 = length(w_update01); 
%pruning, merging, capping 
[w_update01,m_update01,P_update01]= gaus_prune(w_update01,m_update01,P_update01,filter.elim_threshold);    
L_prune01= length(w_update01);
[w_update01,m_update01,P_update01]= gaus_merge_newnew(w_update01,m_update01,P_update01,filter.merge_threshold);   
L_merge01= length(w_update01);
[w_update01,m_update01,P_update01]= gaus_cap(w_update01,m_update01,P_update01,filter.L_max);
L_cap01 = length(w_update01);
L_update01 = L_cap01;

%CT
L_posterior02 = length(w_update02); 
%pruning, merging, capping 
[w_update02,m_update02,P_update02]= gaus_prune(w_update02,m_update02,P_update02,filter.elim_threshold);    
L_prune02= length(w_update02);
[w_update02,m_update02,P_update02]= gaus_merge_newnew(w_update02,m_update02,P_update02,filter.merge_threshold);   
L_merge02= length(w_update02);
[w_update02,m_update02,P_update02]= gaus_cap(w_update02,m_update02,P_update02,filter.L_max);
L_cap02 = length(w_update02);
L_update02 = L_cap02;


%CT 2
L_posterior03 = length(w_update03); 
%pruning, merging, capping 
[w_update03,m_update03,P_update03]= gaus_prune(w_update03,m_update03,P_update03,filter.elim_threshold);    
L_prune03 = length(w_update03);
[w_update03,m_update03,P_update03]= gaus_merge_newnew(w_update03,m_update03,P_update03,filter.merge_threshold);   
L_merge03= length(w_update03);
[w_update03,m_update03,P_update03]= gaus_cap(w_update03,m_update03,P_update03,filter.L_max);
L_cap03 = length(w_update03);
L_update03 = L_cap03;

% 找出最大的高斯分量数量
max_GMs_mix = max([L_update01,L_update02,L_update03,1]);

% 对 CV 模型进行处理
if L_update01 < max_GMs_mix
    % 填充零权重分量
    num_to_add = max_GMs_mix - L_update01;
    w_update01 = [w_update01; zeros(num_to_add, 1)];
    m_update01 = [m_update01  zeros(model.x_dimCV,num_to_add)];
    P_update01 = cat(3,P_update01, zeros(model.x_dimCV,model.x_dimCV,num_to_add));
end

% 对 CT 模型进行处理
if L_update02 < max_GMs_mix
    num_to_add = max_GMs_mix - L_update02;
    w_update02 = [w_update02; zeros(num_to_add, 1)];
    m_update02 = [m_update02 zeros(model.x_dimCT,num_to_add)];
    P_update02 = cat(3,P_update02, zeros(model.x_dimCT,model.x_dimCT,num_to_add));
end

% 对 CT 模型进行处理
if L_update03 < max_GMs_mix
    num_to_add = max_GMs_mix - L_update03;
    w_update03 = [w_update03; zeros(num_to_add, 1)];
    m_update03 = [m_update03 zeros(model.x_dimCT,num_to_add)];
    P_update03 = cat(3,P_update03, zeros(model.x_dimCT,model.x_dimCT,num_to_add));
end
%% CT
% delete omega
m_update002 = [m_update02([1 2 3 4],:)];
P_update002 = P_update02(1:4, 1:4, :);

%% CT 2
% delete omega
m_update003 = [m_update03([1 2 3 4],:)];
P_update003 = P_update03(1:4, 1:4, :);
end