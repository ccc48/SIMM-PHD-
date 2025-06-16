function model= gen_model

%% basic parameters
model.x_dimCV = 4;  %[x vx y vy]
model.x_dimCT = 5;  %[x vx y vy w]
model.z_dim = 2;    %[x y]dimension of observation vector
model.M=3;          %模型个数 1CV 2CT

model.m=30;          
model.n=30;
model.p=2;
% model.Pmn = [0.98,  0.01, 0.01;
%              0.01, 0.98,  0.01; 
%              0.01, 0.01, 0.98];
% model.Pmn = [0.95,  0.025, 0.025;
%              0.025, 0.95,  0.025; 
%              0.025, 0.025, 0.95];
model.Pmn = [0.9,  0.05, 0.05;%切换√
             0.05, 0.9,  0.05; 
             0.05, 0.05, 0.9];
% model.Pmn = [0.84, 0.08, 0.08;%  1
%              0.08, 0.84, 0.08; 
%              0.08, 0.08, 0.84];
% model.Pmn = [0.8, 0.1, 0.1;%切换√
%              0.1, 0.8, 0.1; 
%              0.1, 0.1, 0.8];

% model.Sij = [  1,  0.4,  0.4; 
%              0.4,    1,  0.4;
%              0.4,  0.4,    1];%运动模型转移可能性矩阵

% model.Sij = [  1,  0.3,  0.3; 
%              0.3,    1,  0.3;
%              0.3,  0.3,    1];%运动模型转移可能性矩阵
model.Sij = [  1,  0.25,  0.25; 
             0.25,    1,  0.25;
             0.25,  0.25,    1];%运动模型转移可能性矩阵
% model.Sij = [  1,  0.2,  0.2; 
%              0.2,    1,  0.2;
%              0.2,  0.2,    1];%运动模型转移可能性矩阵
% model.Sij = [  1,  0.1,  0.1; %切换√
%              0.1,    1,  0.1;
%              0.1,  0.1,    1];%运动模型转移可能性矩阵 
% model.Sij = [  1,  0.08,  0.08; % 1      
%              0.08,    1,  0.08;
%              0.08,  0.08,    1];%运动模型转移可能性矩阵 

% model.Sij = [  1,  0.05, 0.05; %切换√
%              0.05,   1,  0.05;
%              0.05, 0.05,   1];%运动模型转移可能性矩阵
% model.Sij = [  1,  0.025, 0.025; 
%              0.025,   1,  0.025;
%              0.025, 0.025,   1];%运动模型转移可能性矩阵

% model.Sij = [  1,  0.01, 0.01; 
%              0.01,   1,  0.01;
%              0.01, 0.01,   1];%运动模型转移可能性矩阵

%% %% dynamical model parameters (CV model)
%%%%% observation model parameters (noisy x/y only)
model.T= 1;                                     %sampling period
model.A0= [ 1 model.T ; 0 1 ];                         %transition matrix                     
model.F_CV = [model.A0 zeros(2,2);zeros(2,2) model.A0];
model.G0= [ (model.T^2)/2; model.T ];
model.G_CV= [ model.G0 zeros(2,1); zeros(2,1) model.G0];
% model.sigma_cv = 3;%√
model.sigma_cv = 5;
% model.sigma_cv = 6;
% model.sigma_cv = 10;
model.Q_CV= (model.sigma_cv)^2* model.G_CV*model.G_CV';   %process noise covariance

model.H_CV = [1 0 0 0;
              0 0 1 0];
% model.D_CV = diag([ 5; 5 ]);%减小测量标准差至5，原为10
model.D_CV = diag([ 6; 6 ]);
% model.D_CV = diag([ 10; 10 ]);
model.R_CV = model.D_CV*model.D_CV';              %observation noise covariance

%% (CT1 model)
model.omega = -pi/45;
% model.omega = pi/180;  % Turn rate, rad/s 转弯率 %0.01745弧度=1°
a = sin(model.T*model.omega)./model.omega;
b = (1-cos(model.T*model.omega))./model.omega;
c = cos(model.T*model.omega);
d = sin(model.T*model.omega);
model.F_CT = [1   a   0 -b  0; 
              0   c   0 -d  0; 
              0   b   1  a  0; 
              0   d   0  c  0;
              0   0   0  0  1];
model.sigma_vel = 6;
% model.sigma_vel = 20;
% model.sigma_vel = 10;
model.sigma_turn = pi/180;  %大
% model.sigma_turn = pi/720;  %0.25°
model.bt = model.sigma_vel*[ (model.T^2)/2; model.T ];
model.G_CT = [ model.bt zeros(2,2); zeros(2,1) model.bt zeros(2,1); zeros(1,2) model.T*model.sigma_turn ];
model.B_CT = eye(3);
model.Q_CT = model.B_CT*model.B_CT';

model.H_CT = [1 0 0 0 0;        %保留，求meas的时候用
              0 0 1 0 0;
              0 0 0 0 1];
% model.D_CT = diag([ 5; 5]);%减小测量标准差至5，原为10
model.D_CT = diag([ 6; 6 ]);
% model.D_CT = diag([ 10; 10]);%减小测量标准差至5，原为10
model.R_CT = model.D_CT*model.D_CT';              %observation noise covariance

%% (CT2 model)
model.omega2 = pi/60;
% model.omega2 = pi/90;
% model.omega2 = pi/180;  % Turn rate, rad/s 转弯率 %0.01745弧度=1°
a2 = sin(model.T*model.omega2)./model.omega2;
b2 = (1-cos(model.T*model.omega2))./model.omega2;
c2 = cos(model.T*model.omega2);
d2 = sin(model.T*model.omega2);
model.F_CT2 = [1   a2   0 -b2  0; 
               0   c2   0 -d2  0; 
               0   b2   1  a2  0; 
               0   d2   0  c2  0;
               0   0    0  0   1];
% model.sigma_vel2 = 4;
% model.sigma_turn2 = pi/360;  
% model.sigma_vel2 = 15;

% model.sigma_vel2 = 5;
model.sigma_vel2 = 7;
model.sigma_turn2 = pi/360;  
% model.sigma_turn2 = pi/720;  %0.25°
model.bt2 = model.sigma_vel2*[ (model.T^2)/2; model.T ];
model.G_CT2 = [ model.bt2 zeros(2,2); zeros(2,1) model.bt2 zeros(2,1); zeros(1,2) model.T*model.sigma_turn2 ];
model.B_CT2 = eye(3);
model.Q_CT2 = model.B_CT*model.B_CT';

model.H_CT2 = [1 0 0 0 0;        %保留，求meas的时候用
               0 0 1 0 0;
               0 0 0 0 1];
% model.D_CT2 = diag([ 5; 5]);%减小测量标准差至5，原为10
model.D_CT2 = diag([ 6; 6 ]);
% model.D_CT2 = diag([ 10; 10]);%减小测量标准差至5，原为10
model.R_CT2 = model.D_CT2*model.D_CT2';              %observation noise covariance

model.G = cell(2, 1); 
model.G{1} = model.G_CV;
model.G{2} = model.G_CT;
model.G{3} = model.G_CT2;

model.Q = cell(2, 1); 
model.Q{1} = model.Q_CV;
model.Q{2} = model.Q_CT;
model.Q{3} = model.Q_CT2;

% %出生死亡时间
% % Define maneuver intervals and models
model.maneuvers = [1, 61, 101; 1, 2, 3];  % Start time and model index (1-CV, 2-CA, 3-CT)
model.tbirth = [1, 61, 101];
model.tdeath = [60,100,150];
% model.maneuvers = [1, 51, 101; 1, 2, 3];  % Start time and model index (1-CV, 2-CA, 3-CT)
% model.tbirth = [1, 51, 101];
% model.tdeath = [50,100,150];

%% survival/death parameters
model.P_S= .99;
model.Q_S= 1-model.P_S;

%% birth parameters (Poisson birth model, multiple Gaussian components)
model.L_birth= 4;                                                     %no. of Gaussian birth terms
model.w_birth= zeros(model.L_birth,1);                                %weights of Gaussian birth terms (per scan) 
                                                                      % [sum gives average rate of target birth per scan]
%CV
model.m_birthCV = zeros(model.x_dimCV,model.L_birth);                      %means of Gaussian birth terms 
model.B_birthCV = zeros(model.x_dimCV,model.x_dimCV,model.L_birth);          %std of Gaussian birth terms
model.P_birthCV = zeros(model.x_dimCV,model.x_dimCV,model.L_birth);          %cov of Gaussian birth terms
%CT
model.m_birthCT = zeros(model.x_dimCT,model.L_birth);                      %means of Gaussian birth terms 
model.B_birthCT = zeros(model.x_dimCT,model.x_dimCT,model.L_birth);          %std of Gaussian birth terms
model.P_birthCT = zeros(model.x_dimCT,model.x_dimCT,model.L_birth);          %cov of Gaussian birth terms

%CT 2
model.m_birthCT2 = zeros(model.x_dimCT,model.L_birth);                      %means of Gaussian birth terms 
model.B_birthCT2 = zeros(model.x_dimCT,model.x_dimCT,model.L_birth);          %std of Gaussian birth terms
model.P_birthCT2 = zeros(model.x_dimCT,model.x_dimCT,model.L_birth);          %cov of Gaussian birth terms

%% %%birth term 1 ~%%%%
model.w_birth(1)= 3/100;                                              
%% term 1 ~ CV
% model.m_birthCV(:,1)= [  1300; 25;  1200; 0];
% model.m_birthCV(:,1)= [  1603; -20;  1140; 34];%t=60
model.m_birthCV(:,1)= [  1583; -20;  1175; 34];%t=61
% model.m_birthCV(:,1)= [  1866; 0;  1702; 0];%t=51

model.B_birthCV(:,:,1)= diag([ 10; 1; 10; 1 ]);
% model.B_birthCV(:,:,1)= diag([ 5; 1; 5; 1 ]);
% model.B_birthCV(:,:,1)= diag([ 10; 10; 10; 10 ]);

model.P_birthCV(:,:,1)= model.B_birthCV(:,:,1)*model.B_birthCV(:,:,1)';

%%  term 1 ~ CT
% model.m_birthCT(:,1)= [  1300; 25;  1200; 0; -pi/45 ];
% model.m_birthCT(:,1)= [2100; 40;  2000; 0; -pi/45 ];%xstart
model.m_birthCT(:,1)= [2140; 40;  1998; -2.79; -pi/45 ];%t=1

model.B_birthCT(:,:,1)= diag([ 10; 5; 10; 5; 3*(pi/180) ]);
% model.B_birthCT(:,:,1)= diag([ 10; 10; 10; 10; 6*(pi/720) ]);
model.P_birthCT(:,:,1)= model.B_birthCT(:,:,1)*model.B_birthCT(:,:,1)';

%%  term 1 ~ CT 2
% model.m_birthCT2(:,1)= [  1300; 25;  1200; 0; -pi/45 ];
% model.m_birthCT2(:,1)= [803; -20;  2526; 34; -pi/45 ];%t=100
model.m_birthCT2(:,1)= [783; -21;  2560; 34; -pi/60 ];%t=101
% model.m_birthCT2(:,1)= [-13; 0;  2385; 0; pi/60 ];%t=101

model.B_birthCT2(:,:,1)= diag([ 10; 5; 10; 5; 3*(pi/180) ]);
% model.B_birthCT2(:,:,1)= diag([ 10; 10; 10; 10; 6*(pi/720) ]);
model.P_birthCT2(:,:,1)= model.B_birthCT(:,:,1)*model.B_birthCT(:,:,1)';


%% %%birth term 2 ~%%%%
model.w_birth(2)= 3/100;                                              
%% term 2 ~ CV
% model.m_birthCV(:,2)= [ 400; 0; -600; 0 ];
model.m_birthCV(:,2)= [ 2000; 0; -3000; 0 ];

model.B_birthCV(:,:,2)= diag([ 10; 1; 10; 1 ]);
% model.B_birthCV(:,:,2)= diag([ 5; 1; 5; 1 ]);
% model.B_birthCV(:,:,2)= diag([ 10; 10; 10; 10 ]);
model.P_birthCV(:,:,2)= model.B_birthCV(:,:,2)*model.B_birthCV(:,:,2)';

%% term 2 ~ CT
% model.m_birthCT(:,2)= [ 400; 0; -600; 0; 0 ];
model.m_birthCT(:,2)= [ 2000; 0; -3000; 0; 0 ];

model.B_birthCT(:,:,2)= diag([ 10; 5; 10; 5; 3*(pi/180) ]);
% model.B_birthCT(:,:,2)= diag([ 10; 10; 10; 10; 6*(pi/720) ]);

model.P_birthCT(:,:,2)= model.B_birthCT(:,:,2)*model.B_birthCT(:,:,2)';

%% term 2 ~ CT 2
% model.m_birthCT2(:,2)= [ 400; 0; -600; 0; 0 ];
% model.m_birthCT2(:,2)= [ 2000; 0; -3000; 0; 0 ];
model.m_birthCT2(:,2)= [783; -21;  2560; 34; -pi/60 ];%t=101


model.B_birthCT2(:,:,2)= diag([ 10; 5; 10; 5; 3*(pi/180) ]);
% model.B_birthCT2(:,:,2)= diag([ 10; 10; 10; 10; 6*(pi/720) ]);

model.P_birthCT2(:,:,2)= model.B_birthCT(:,:,2)*model.B_birthCT(:,:,2)';

%% %%birth term 3 ~%%%%
model.w_birth(3)= 3/100;                                              
%% term 3 ~CV
% model.m_birthCV(:,3)= [ -800; 0; -200; 0 ];
model.m_birthCV(:,3)= [ -4000; 0; -1000; 0 ];

model.B_birthCV(:,:,3)= diag([ 10; 1; 10; 1 ]);
% model.B_birthCV(:,:,3)= diag([ 5; 1; 5; 1 ]);
% model.B_birthCV(:,:,3)= diag([ 10; 10; 10; 10 ]);

model.P_birthCV(:,:,3)= model.B_birthCV(:,:,3)*model.B_birthCV(:,:,3)';

%% term 3 ~ CT
% model.m_birthCT(:,3)= [ -800; 0; -200; 0; 0 ]; 
model.m_birthCT(:,3)= [ -4000; 0; -1000; 0; 0 ]; 

model.B_birthCT(:,:,3)= diag([10; 5; 10; 5; 3*(pi/180) ]);
% model.B_birthCT(:,:,3)= diag([ 10; 10; 10; 10; 6*(pi/720) ]);

model.P_birthCT(:,:,3)= model.B_birthCT(:,:,3)*model.B_birthCT(:,:,3)';

%% term 3 ~ CT 2
% model.m_birthCT2(:,3)= [ -800; 0; -200; 0; 0 ]; 
model.m_birthCT2(:,3)= [ -4000; 0; -1000; 0; 0 ]; 

model.B_birthCT2(:,:,3)= diag([10; 5; 10; 5; 3*(pi/180) ]);
% model.B_birthCT2(:,:,3)= diag([ 10; 10; 10; 10; 6*(pi/720) ]);

model.P_birthCT2(:,:,3)= model.B_birthCT(:,:,3)*model.B_birthCT(:,:,3)';

%% %%birth term 4 ~%%%%
model.w_birth(4)= 3/100;                                              
%% term 4 ~ CV
% model.m_birthCV(:,4)= [ -200; 0; 800; 0 ];
model.m_birthCV(:,4)= [ -1000; 0; 4000; 0 ];

model.B_birthCV(:,:,4)= diag([ 10; 1; 10; 1 ]);
% model.B_birthCV(:,:,4)= diag([ 5; 1; 5; 1 ]);
% model.B_birthCV(:,:,4)= diag([ 10; 10; 10; 10 ]);

model.P_birthCV(:,:,4)= model.B_birthCV(:,:,4)*model.B_birthCV(:,:,4)';
%% term 4 ~ CT
% model.m_birthCT(:,4)= [ -200; 0; 800; 0; 0 ];
model.m_birthCT(:,4)= [ -1000; 0; 4000; 0; 0 ];
% model.m_birthCT(:,4)= [ 2669; 0; 1486; 0; 0 ]; %！

model.B_birthCT(:,:,4)= diag([ 10; 5; 10; 5; 3*(pi/180) ]);
% model.B_birthCT(:,:,4)= diag([ 10; 10; 10; 10; 6*(pi/720) ]);
model.P_birthCT(:,:,4)= model.B_birthCT(:,:,4)*model.B_birthCT(:,:,4)';

%% term 4 ~ CT 2
% model.m_birthCT2(:,4)= [ -200; 0; 800; 0; 0 ];
model.m_birthCT2(:,4)= [ -1000; 0; 4000; 0; 0 ];

model.B_birthCT2(:,:,4)= diag([ 10; 5; 10; 5; 3*(pi/180) ]);
% model.B_birthCT2(:,:,4)= diag([ 10; 10; 10; 10; 6*(pi/720) ]);
model.P_birthCT2(:,:,4)= model.B_birthCT(:,:,4)*model.B_birthCT(:,:,4)';



% detection parameters
% model.P_D= 0.999;   %probability of detection in measurements
model.P_D= 0.99;   %probability of detection in measurements
% model.P_D= 0.98;   %probability of detection in measurements 1

% model.P_D= 0.95;   %probability of detection in measurements
% model.P_D= 0.90;   %probability of detection in measurements pd越小，SIMM的性能越差
model.Q_D= 1-model.P_D; %probability of missed detection in measurements

% clutter parameters
model.lambda_c= 60; 
% model.lambda_c= 20; 
% model.lambda_c= 10; % 1 减少杂波平均数量，改为10（原为20） %泊松分布的均值，表示每个扫描周期内的平均杂波数%poisson average rate of uniform clutter (per scan)
% model.range_c= [ -1000 2000; 0 3000 ];%二维区域的上下界      
model.range_c= [ -1500 3500; 0 4000 ];%二维区域的上下界      

model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
%均匀杂波的密度，计算方式是1除以区域面积，因为每个维度的范围差相乘得到面积，所以概率密度就是1/面积。
%prod(model.range_c(:,2)-model.range_c(:,1))，这在二维得到的是长乘宽，

end