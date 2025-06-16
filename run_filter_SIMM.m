function estS = run_filter_SIMM(model,meas,truth)
%=== Setup

%% 3个模型各自的output variables
estS.W= cell(meas.K,model.M);
estS.X= cell(meas.K,model.M);       
estS.P= cell(meas.K,model.M);
estS.N= zeros(meas.K,model.M);
estS.L= cell(meas.K,model.M);      %轨迹数？
estS.Lk= cell(meas.K,model.M);     %似然
%% %交互输出的
estS.IMMW = cell(meas.K,1);
estS.IMMX = cell(meas.K,1);    
estS.IMMP = cell(meas.K,1);
estS.IMMN = zeros(meas.K,1);
estS.IMML = cell(meas.K,1);
estS.IMMLk = cell(meas.K,1);

%filter parameters
filter.L_max= 100;                  %limit on number of Gaussians
filter.elim_threshold= 1e-4;        %pruning threshold
filter.merge_threshold= 4;          %merging threshold
filter.P_G= 0.999;                               %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0
filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

estS.filter= filter;

%=== Filtering 

%initial prior
w_update(1,:)= eps;
m_update(:,1)= [0.1; 0;  
                0.1; 0];
P_update(:,:,1) = diag([10 1 10 1]).^2;
L_update = 1;

       
miu=[1; 1; 1];
estS.miu = [];
%% recursive filtering
for k=1:meas.K
    %强制切换
    % if k==61
    %     miu = [0.8; 0.1; 0.1];
    % elseif k==101
    %     miu = [0.1; 0.8; 0.1];
    % end
    
    %% 输入交互 %% input interaction
    mu = zeros(model.M,model.M);    % 混合概率矩阵
    c = zeros(model.M,1);     % 归一化常数C_j
    for j = 1:model.M
        c(j,:) = max(model.Sij(:,j).*miu);% 计算归一化常数 (使用max而非IMM中的sum)
    end
    for i = 1:model.M
        for j = 1:model.M
            mu(i,j) = model.Sij(i,j) * miu(i) / c(j);% 计算混合概率
        end
    end
    if k==1
        % 初始时刻使用相同的初始值
        w_1M = w_update;  w_2M = w_update;  w_3M = w_update;  % 获取前一时刻各模型的状态
        m_1M = m_update;  m_2M = m_update;  m_3M = m_update;  
        P_1M = P_update;  P_2M = P_update;  P_3M = P_update;  
    else
        % 非初始时刻使用前一时刻估计值
        w_1M = estS.W{k-1,1};  w_2M = estS.W{k-1,2};  w_3M = estS.W{k-1,3};  
        m_1M = estS.X{k-1,1};  m_2M = estS.X{k-1,2};  m_3M = estS.X{k-1,3}; 
        P_1M = estS.P{k-1,1};  P_2M = estS.P{k-1,2};  P_3M = estS.P{k-1,3}; 
    end
    WW = {w_1M  w_2M  w_3M };
    MM = {m_1M  m_2M  m_3M };
    PP = {P_1M  P_2M  P_3M };
    % 计算交互输入 - 加权混合所有模型的状态
    for j = 1:model.M
        [~,lmax] = max(mu(:,j));  %
        w0{j} = WW{lmax};        %权重、状态
        m0{j} = MM{lmax};
    end
    for j = 1:model.M    
        [~,l] = max(mu(:,j));  
        P0{j} = PP{j} + (MM{j} - m0{l})*(MM{j} - m0{l})';%协方差
    end
    L_update1 = length(w0{1});
    L0 = {L_update1;L_update1;L_update1};
    
    % 滤波并计算似然
    [w_update01,m_update01,P_update01,L_update01,LK01,...
     w_update02,m_update02,P_update02,L_update02,LK02,...
     w_update03,m_update03,P_update03,L_update03,LK03] = CVCTfilters(k,model,meas,estS,filter,w0,m0,P0,L0);

    estS.W{k,1} = w_update01;estS.W{k,2} = w_update02;estS.W{k,3} = w_update03;
    estS.X{k,1} = m_update01;estS.X{k,2} = m_update02;estS.X{k,3} = m_update03;
    estS.P{k,1} = P_update01;estS.P{k,2} = P_update02;estS.P{k,3} = P_update03;
    estS.L{k,1} = L_update01;estS.L{k,2} = L_update02;estS.L{k,3} = L_update03;
    estS.Lk{k,1} = LK01;     estS.Lk{k,2} = LK02;     estS.Lk{k,3} = LK03;

    %% ---更新模型概率
    estS.miu = [estS.miu miu];
    Lk = [estS.Lk{k,1}  estS.Lk{k,2}  estS.Lk{k,3}];
    %找到Lk的最小值
    if all(Lk(:) == 0)     %全零矩阵
        Lk_min = eps;  % 指定默认值（如 eps）
    else
        Lk_min = min(Lk(Lk ~= 0));  %找出Lk矩阵中非零元素的最小值
    end
    for mll =1:model.M       
        if Lk(mll) == 0
            Lk(mll) = Lk_min*eps;
        end
    end
    c_norm = max(Lk .* c');  %归一化常数

    for mll =1:model.M
        miu(mll,:) = Lk(:,mll) * c(mll,:)./c_norm;
    end
    % miu = miu/max(miu);

    %% 输出交互
    [w_update,m_update,P_update,L_update] = SIMM_output_interaction(k,model,estS,miu);
          
    % --- mixture management
    L_posterior= length(w_update);
    
    %pruning, merging, capping
    [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,filter.elim_threshold);    L_prune= length(w_update);
    [w_update,m_update,P_update]= gaus_merge_newnew(w_update,m_update,P_update,filter.merge_threshold);   L_merge= length(w_update);
    [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,filter.L_max);
    L_cap  = length(w_update);
    
    L_update= L_cap;
    
    %--- state extraction
    idx= find(w_update > 0.5 );
    for j=1:length(idx)
        repeat_num_targets= round(w_update(idx(j)));
        estS.IMMW{k}= [ estS.IMMW{k} repmat(w_update(idx(j),:),[1,repeat_num_targets]) ];
        estS.IMMX{k}= [ estS.IMMX{k}  repmat(m_update(:,idx(j)),[1,repeat_num_targets]) ];
        estS.IMMN(k)= estS.IMMN(k)+repeat_num_targets;
        estS.IMML{k}= [];
    end

end

for k=1:meas.K
    % 检查est.IMMX{k}是否为空
    if ~isempty(estS.IMMX{k})  %非空
        estS.error_X{k} = estS.IMMX{k} - truth.X{k};
    else
        estS.IMMX{k} = [0;0;0;0];  %补零
        estS.error_X{k} = estS.IMMX{k} - truth.X{k};%看误差
        estS.IMMX{k} = [];%返回空矩阵，不然绘图的时候位置不对
    end
end

end