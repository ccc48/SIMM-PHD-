function est = run_filter_IMM(model,meas,truth)
%=== Setup

%% 模型各自的 output variables
est.W= cell(meas.K,model.M);
est.X= cell(meas.K,model.M);        
est.P= cell(meas.K,model.M);
est.N= zeros(meas.K,model.M);      % 势（数目）
est.L= cell(meas.K,model.M);       % 轨迹数
est.Lk= cell(meas.K,model.M);      % 似然
%% 交互输出的output variables
est.IMMW = cell(meas.K,1);
est.IMMX = cell(meas.K,1);   
est.IMMP = cell(meas.K,1);
est.IMMN = zeros(meas.K,1);
est.IMML = cell(meas.K,1);
est.IMMLk = cell(meas.K,1);

%%滤波参数filter parameters
filter.L_max= 100;                  % limit on number of Gaussians
filter.elim_threshold= 1e-4;        % pruning threshold%值越大，高斯分量越容易被剪枝
                                    % 1e-3太大了，有可能所有高斯分量被剪枝
filter.merge_threshold= 4;          % merging threshold值越大越容易合并
                                    % 任意两个高斯分量的马氏距离小于merge时被合并

filter.P_G= 0.999;                               % gate size in percentage
% filter.P_G= 0.9999999999999999;                               % gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   % inv chi^2 dn gamma value
filter.gate_flag= 1;                             % gating on or off 1/0
filter.run_flag= 'disp';            % 'disp' or 'silence' for on the fly output
est.filter= filter;


%% initial prior
w_update(1,:)= eps;
m_update(:,1)= [0.1; 0; 
                0.1; 0];
P_update(:,:,1) = diag([10 1 10 1]).^2;
L_update = 1;

miu=[1/3; 1/3; 1/3];
est.miu = [];    % 模型概率
%% === Filtering 

%% recursive filtering
for k=1:meas.K
    %强制切换
    % if k==61
    %     miu = [0.8; 0.1; 0.1];
    % elseif k==101
    %     miu = [0.1; 0.8; 0.1];
    % end
    %% 输入交互 %% input interaction
    mu = zeros(model.M,model.M);    %混合概率矩阵
    c = zeros(model.M,1);     %归一化常数C_j
    for j = 1:model.M
        c(j,:) =  sum( model.Pmn(:,j).*miu );
    end
    for i = 1:model.M
        for j = 1:model.M
            mu(i,j) = model.Pmn(i,j) *miu(i)./c(j);
        end
    end
    if k==1
        w_1M = w_update;  w_2M = w_update;  w_3M = w_update;  
        m_1M = m_update;  m_2M = m_update;  m_3M = m_update;  
        P_1M = P_update;  P_2M = P_update;  P_3M = P_update; 
        N_1M = L_update;  N_2M = L_update;  N_3M = L_update; 
    else
        w_1M = est.W{k-1,1};  w_2M = est.W{k-1,2};  w_3M = est.W{k-1,3}; 
        m_1M = est.X{k-1,1};  m_2M = est.X{k-1,2};  m_3M = est.X{k-1,3};  
        P_1M = est.P{k-1,1};  P_2M = est.P{k-1,2};  P_3M = est.P{k-1,3};
    end
    % 计算每个CV分量与每个CT分量和每个CA分量的交互输入
    for j = 1:model.M
        w0{j} = w_1M.*mu(1,j) + w_2M.*mu(2,j) + w_3M.*mu(3,j) ;
        m0{j} = m_1M.*mu(1,j) + m_2M.*mu(2,j) + m_3M.*mu(3,j) ;
    end
    for j = 1:model.M
        P0{j} = mu(1,j) * (P_1M + (m_1M - m0{j}) * (m_1M - m0{j})')... 
               +mu(2,j) * (P_2M + (m_2M - m0{j}) * (m_2M - m0{j})')...
               +mu(3,j) * (P_3M + (m_3M - m0{j}) * (m_3M - m0{j})');
            
    end
    L_update1 = length(w0{1});
    L0 = {L_update1;L_update1;L_update1};
    
    [w_update01,m_update01,P_update01,L_update01,LK01,...
     w_update02,m_update02,P_update02,L_update02,LK02,...
     w_update03,m_update03,P_update03,L_update03,LK03] = CVCTfilters(k,model,meas,est,filter,w0,m0,P0,L0);

    est.W{k,1} = w_update01;est.W{k,2} = w_update02;est.W{k,3} = w_update03;
    est.X{k,1} = m_update01;est.X{k,2} = m_update02;est.X{k,3} = m_update03;
    est.P{k,1} = P_update01;est.P{k,2} = P_update02;est.P{k,3} = P_update03;
    est.L{k,1} = L_update01;est.L{k,2} = L_update02;est.L{k,3} = L_update03;
    est.Lk{k,1} = LK01;     est.Lk{k,2} = LK02;     est.Lk{k,3} = LK03;
    
    %% ---更新模型概率 model probability update
    est.miu = [est.miu miu];
    Lk = [est.Lk{k,1}  est.Lk{k,2}  est.Lk{k,3} ];
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
    c_norm = Lk * c;  %归一化常数  %c_norm和c意义完全不一样，注意

    for mll =1:model.M
        miu(mll,:) = Lk(:,mll) * c(mll,:)./c_norm;
    end

    %% 输出交互
    [w_update,m_update,P_update,L_update] = IMM_output_interaction(k,model,est,miu);


            
    %% --- mixture management
    L_posterior = length(w_update);
    
    %% pruning, merging, capping
    [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,filter.elim_threshold);    L_prune= length(w_update);
    [w_update,m_update,P_update]= gaus_merge_newnew(w_update,m_update,P_update,filter.merge_threshold);   L_merge= length(w_update);
    [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,filter.L_max);
    L_cap  = length(w_update);
    
    L_update= L_cap;
    
    %% --- state extraction
    idx= find(w_update > 0.5 );
    for j=1:length(idx)
        repeat_num_targets= round(w_update(idx(j)));
        est.IMMW{k}= [ est.IMMW{k} repmat(w_update(idx(j),:),[1,repeat_num_targets]) ];
        est.IMMX{k}= [ est.IMMX{k}  repmat(m_update(:,idx(j)),[1,repeat_num_targets]) ];
        est.IMMN(k)= est.IMMN(k)+repeat_num_targets;
        est.IMML{k}= [];
    end
   
end

%% RMSE
for k=1:meas.K
    % 检查est.IMMX{k}是否为空
    if ~isempty(est.IMMX{k})  %非空
        est.error_X{k} = est.IMMX{k} - truth.X{k};
    else
        est.IMMX{k} = [0;0;0;0];  %补零
        est.error_X{k} = est.IMMX{k} - truth.X{k};
        est.IMMX{k} = []; %返回空矩阵，不然绘图的时候位置不对
    end
end

end

            