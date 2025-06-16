function X= gen_newstate_fn(k,model,Xd,G,Q,V)

% Determining the correct noise parameter sigma_V based on the matrix F
%% CV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isequal(G, model.G_CV)               
        sigma = model.sigma_cv;
        G = model.G_CV;
         % Generate noise based on the model   
        if ~isnumeric(V)
            if strcmp(V,'noise')
                V= sigma*G*randn(size(G,2),size(Xd,2));
            elseif strcmp(V,'noiseless')
                V= zeros(size(G,1),size(Xd,2));
            end
        end

        % Generate new state
        if isempty(Xd)
            X= [];
        else
            X= model.F_CV*Xd+ V;
        end

%% CT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif isequal(G, model.G_CT)
        %nonlinear state space equation (CT model)
        % 4--> 5维
        X_1i = [Xd([1 2 3 4],:); model.omega]; % 加入OMEGA，构成5*n矩阵
        
        if ~isnumeric(V)
            if strcmp(V,'noise')
                V= model.B_CT*randn(size(model.B_CT,2),size(Xd,2));
            elseif strcmp(V,'noiseless')
                V= zeros(size(model.B_CT,1),size(Xd,2));
            end
        end

        if isempty(Xd)
            X= [];
        else %modify below here for user specified transition model
            % X_2i = model.F_CT*X_1i+ V;          %直接用状态转移矩阵 %第一种方法

            X_2i = zeros(size(X_1i));                   %以下为第二种方法
            %-- short hand
            L= size(X_1i,2);
            T= model.T; 
            omega= model.omega;
            tol= 1e-10;
            %-- pre calcs
            sin_omega_T= sin(omega*T);
            cos_omega_T= cos(omega*T);
            a= T*ones(1,L); b= zeros(1,L);
            idx= find( abs(omega) > tol );
            a(idx)= sin_omega_T(idx)./omega(idx);
            b(idx)= (1-cos_omega_T(idx))./omega(idx);
            %--  x/y pos/vel
            X_2i(1,:)= X_1i(1,:)+ a.*X_1i(2,:)- b.*X_1i(4,:);
            X_2i(2,:)= cos_omega_T.*X_1i(2,:)- sin_omega_T.*X_1i(4,:);
            X_2i(3,:)= b.*X_1i(2,:) + X_1i(3,:)+ a.*X_1i(4,:);
            X_2i(4,:)= sin_omega_T.*X_1i(2,:)+ cos_omega_T.*X_1i(4,:);
            %-- turn rate
            X_2i(5,:)= X_1i(5,:);
            %-- add scaled noise 
            X_2i= X_2i+ model.G_CT*V;

            % delete omega
            X = [X_2i([1 2 3 4],:)];
            %第三种：找到圆心位置、初始相位、初始速度
            %% 推导CT运动参数
            % v = sqrt(X_1i(2,:)^2 + X_1i(4,:)^2);    % 速度大小
            % r = v / abs(model.omega);         % 圆周运动半径
            % xc = X_1i(1,:) - vy0/model.omega;        % 圆心x坐标
            % yc = X_1i(3,:) + vx0/model.omega;        % 圆心y坐标
            % phi0 = atan2(X_1i(3,:) - yc, X_1i(1,:) - xc); % 初始相位
            % 
            % % CT运动轨迹
            % X(1,:)  = xc + r * cos(model.omega*model.T + phi0);
            % X(3,:)  = yc + r * sin(model.omega*model.T + phi0);



        end
        %% CT2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif isequal(G, model.G_CT2)
        %nonlinear state space equation (CT model)
        % 4--> 5维
        X_1i = [Xd([1 2 3 4],:); model.omega2]; % 加入OMEGA，构成5*n矩阵
        
        if ~isnumeric(V)
            if strcmp(V,'noise')
                V= model.B_CT2*randn(size(model.B_CT2,2),size(Xd,2));
            elseif strcmp(V,'noiseless')
                V= zeros(size(model.B_CT2,1),size(Xd,2));
            end
        end

        if isempty(Xd)
            X= [];
        else %modify below here for user specified transition model
            % X_2i = model.F_CT*X_1i+ V;          %直接用状态转移矩阵 %第一种方法

            X_2i = zeros(size(X_1i));                   %以下为第二种方法
            %-- short hand
            L= size(X_1i,2);
            T= model.T; 
            omega= model.omega2;
            tol= 1e-10;
            %-- pre calcs
            sin_omega_T= sin(omega*T);
            cos_omega_T= cos(omega*T);
            a= T*ones(1,L); b= zeros(1,L);
            idx= find( abs(omega) > tol );
            a(idx)= sin_omega_T(idx)./omega(idx);
            b(idx)= (1-cos_omega_T(idx))./omega(idx);
            %--  x/y pos/vel
            X_2i(1,:)= X_1i(1,:)+ a.*X_1i(2,:)- b.*X_1i(4,:);
            X_2i(2,:)= cos_omega_T.*X_1i(2,:)- sin_omega_T.*X_1i(4,:);
            X_2i(3,:)= b.*X_1i(2,:) + X_1i(3,:)+ a.*X_1i(4,:);
            X_2i(4,:)= sin_omega_T.*X_1i(2,:)+ cos_omega_T.*X_1i(4,:);
            %-- turn rate
            X_2i(5,:)= X_1i(5,:);
            %-- add scaled noise 
            X_2i= X_2i+ model.G_CT2*V;

            % delete omega
            X = [X_2i([1 2 3 4],:)];
    

        end
    else
        error('Unrecognized state transition matrix F.');
    end
end