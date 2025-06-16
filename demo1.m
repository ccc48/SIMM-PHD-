clc
clear
tic
MC=100;                   %序列蒙特卡洛次数
ospa_c= 30;
ospa_p= 1;

model =  gen_model;              %model参数 
truth =  gen_truth(model);       %和truth每一次都是一样的，不用变

%% ————————————————————————————————————————————————————————————————
[X_track,k_birth,k_death]= extract_tracks(truth.X,truth.track_list,truth.total_tracks);
%% plot ground truths  
limit= [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) ];
figure; truths= gcf; hold on;
for i=1:truth.total_tracks
    Pt= X_track(:,k_birth(i):1:k_death(i),i); Pt=Pt([1 3],:);%x,y
    plot( Pt(1,[model.tbirth(1):model.tdeath(1)]),Pt(2,[model.tbirth(1):model.tdeath(1)]),'b-'); 
    plot( Pt(1,[model.tbirth(2):model.tdeath(2)]),Pt(2,[model.tbirth(2):model.tdeath(2)]),'r-'); 
    plot( Pt(1,[model.tbirth(3):model.tdeath(3)]),Pt(2,[model.tbirth(3):model.tdeath(3)]),'k-'); 
    
    
    plot( Pt(1,model.tbirth(1)), Pt(2,model.tbirth(1)), 'bo','MarkerSize',6);
    plot( Pt(1,model.tbirth(2)), Pt(2,model.tbirth(2)), 'ro','MarkerSize',6);
    plot( Pt(1,model.tbirth(3)), Pt(2,model.tbirth(3)), 'ko','MarkerSize',6);
    plot( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)), 'k^','MarkerSize',6);
end
axis equal; axis(limit); title('Ground Truths');
xlabel('x-coordinate (m)'); ylabel('y-coordinate (m)');
%% ————————————————————————————————————————————————————————————————

%% 创建进度条窗口 
h = waitbar(0, '正在处理，请稍候...'); 
for MCell =1:MC
    meas= gen_meas(model,truth);    %每次测量随机

    %% IMM-PHD滤波
    est = run_filter_IMM(model,meas,truth);   
    est_MC_X (MCell,:)= est.IMMX;      %滤波值
    ospa_vals= zeros(truth.K,3);
    for k=1:meas.K
        [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est.IMMX{k,:},[1 3]),ospa_c,ospa_p);
    end
    ospa_vals1(:,MCell)=ospa_vals(:,1);
    ospa_vals2(:,MCell)=ospa_vals(:,2);
    ospa_vals3(:,MCell)=ospa_vals(:,3);
    est_MC_N (:,MCell)= est.IMMN;   %势估计
    est_MC_miu {MCell,:} = est.miu;           %模型概率
    est_MC_error (MCell,:)= est.error_X;      %误差值
    

    %% SIMM-PHD滤波
    estS = run_filter_SIMM(model,meas,truth);  
    est_MC_SIMMX (MCell,:)= estS.IMMX;      %滤波值
    ospa_vals_SIMM= zeros(truth.K,3);
    for k=1:meas.K
        [ospa_vals_SIMM(k,1), ospa_vals_SIMM(k,2), ospa_vals_SIMM(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(estS.IMMX{k,:},[1 3]),ospa_c,ospa_p);
    end
    ospa_vals1_SIMM(:,MCell)=ospa_vals_SIMM(:,1);
    ospa_vals2_SIMM(:,MCell)=ospa_vals_SIMM(:,2);
    ospa_vals3_SIMM(:,MCell)=ospa_vals_SIMM(:,3);
    est_MC_N_SIMM (:,MCell)= estS.IMMN;       %势估计
    est_MC_miu_SIMM {MCell,:} = estS.miu;         %模型可能性
    est_MC_error_SIMM (MCell,:)= estS.error_X;      %误差值
    


    %% 更新进度条
    waitbar(MCell/MC, h, sprintf('已完成 %d%%', round(MCell/MC*100)));
end
%% ————————————————————————————————————————————————————————————————
%% 初始化一个空的IMM-滤波值
est.X_aver=cell(truth.K,1);         
[p,q]=size(est_MC_X);
for i=1:q
    xx = [];
    for  l=1:MC
    xx1 = est_MC_X (l,i);
    xx = [xx , cell2mat(xx1)];
    end
    aver_xx=sum(xx,2)/MC;         %貌似是求MC次滤波值的均值
    est.X_aver{i} =[est.X_aver{i} aver_xx];

end

%%
estS.IMM_X_aver=cell(truth.K,1);         %初始化一个空的IMM-ber滤波值
[p_SIMM,q_SIMM]=size(est_MC_SIMMX);
for i=1:q_SIMM
    xx = [];
    for  l=1:MC
    xx1 = est_MC_SIMMX (l,i);
    xx = [xx , cell2mat(xx1)];
    end
    aver_xx_SIMM=sum(xx,2)/MC;         %貌似是求MC次滤波值的均值
    estS.IMM_X_aver{i} =[estS.IMM_X_aver{i} aver_xx_SIMM];

end

%% 绘制测量值、真实值、估计值
%%plot tracks and measurements in x/y
figure; tracking= gcf; hold on;
subplot(211); box on; 

%plot x measurement
for k=1:meas.K
    if ~isempty(meas.Z{k})
        hlined= line(k*ones(size(meas.Z{k},2),1),meas.Z{k}(1,:),'LineStyle','none','Marker','x','Markersize',5,'Color',0.7*ones(1,3));
    end   
end

% plot x track
for i=1:truth.total_tracks
    Px= X_track(:,k_birth(i):1:k_death(i),i); Px=Px([1 3],:);
    hline1= line(k_birth(i):1:k_death(i),Px(1,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
end

% plot x estimate(IMM
for k=1:meas.K
    if ~isempty(est.X_aver{k})
        P= est.X_aver{k}([1 3],:);
        hline2= line(k*ones(size(est.X_aver{k},2),1),P(1,:),'LineStyle','none','Marker','.','Markersize',9,'Color','k');
    end
end

% plot x estimate(SIMM
for k=1:meas.K
    if ~isempty(estS.IMM_X_aver{k})
        P= estS.IMM_X_aver{k}([1 3],:);
        hline3= line(k*ones(size(estS.IMM_X_aver{k},2),1),P(1,:),'LineStyle','none','Marker','.','Markersize',8,'Color','b');
    end
end

% #y#
subplot(212); box on;
%plot y measurement  
for k=1:meas.K
    if ~isempty(meas.Z{k})
        yhlined= line(k*ones(size(meas.Z{k},2),1),meas.Z{k}(2,:),'LineStyle','none','Marker','x','Markersize',5,'Color',0.7*ones(1,3));
    end
end

% plot y track
for i=1:truth.total_tracks
        Py= X_track(:,k_birth(i):1:k_death(i),i); Py=Py([1 3],:);
        yhline1= line(k_birth(i):1:k_death(i),Py(2,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
end

%plot y estimate(IMM
for k=1:meas.K
    if ~isempty(est.X_aver{k})
        P= est.X_aver{k}([1 3],:);
        yhline2= line(k*ones(size(est.X_aver{k},2),1),P(2,:),'LineStyle','none','Marker','.','Markersize',9,'Color','k');
    end
end

% %plot y estimate(SIMM
for k=1:meas.K
    if ~isempty(estS.IMM_X_aver{k})
        P= estS.IMM_X_aver{k}([1 3],:);
        yhline3= line(k*ones(size(estS.IMM_X_aver{k},2),1),P(2,:),'LineStyle','none','Marker','.','Markersize',8,'Color','b');
    end
end

subplot(211); xlabel('Time'); ylabel('x-coordinate (m)');
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(1,:));
legend([hline3 hline2 hline1 hlined],'SIMM Estimates','IMM Estimates','True tracks','Measurements');

subplot(212); xlabel('Time'); ylabel('y-coordinate (m)');
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(2,:));
% legend([yhline2 yhline1 yhlined],'Estimates          ','True tracks','Measurements');


%% OSPA --- Dist Loc Card 
ospa_1=sum(ospa_vals1,2)/MC;   %取MC次的均值
ospa_2=sum(ospa_vals2,2)/MC;
ospa_3=sum(ospa_vals3,2)/MC;

ospa_1SIMM=sum(ospa_vals1_SIMM,2)/MC;   %取MC次的均值
ospa_2SIMM=sum(ospa_vals2_SIMM,2)/MC;
ospa_3SIMM=sum(ospa_vals3_SIMM,2)/MC;

%% 势估计
est.N_aver=zeros(truth.K,1);         %初始化
[p,q]=size(est_MC_N);
for i=1:p
    nn = [];
    for  l=1:MC
    nn1 = est_MC_N (i,l);
    nn = [nn , nn1];
    end
    aver_nn=sum(nn,2)/MC;         %貌似是求MC次滤波值的均值
    est.N_aver(i,:) = round(aver_nn);
end

estS.N_aver=zeros(truth.K,1);         %初始化
[p_SIMM,q_SIMM]=size(est_MC_N_SIMM);
for i=1:p_SIMM
    nn = [];
    for  l=1:MC
    nn1 = est_MC_N_SIMM (i,l);
    nn = [nn , nn1];
    end
    aver_nn_SIMM=sum(nn,2)/MC;         %貌似是求MC次滤波值的均值
    estS.N_aver(i,:) = round(aver_nn_SIMM);

end

%% miu  IMM模型概率切换
est.miu_aver = [];         %初始化一个空的IMM-ber滤波值
[p,q]=size(est_MC_miu{MC,:});
for i=1:p        %3种模型
    xx = [];     %
    for  l=1:MC  %MC次循环
    xx1 = est_MC_miu {l};
    xx = [xx ; xx1(i,:)];
    end
    aver_xx=sum(xx,1)/MC;         %貌似是求MC次滤波值的均值
    est.miu_aver(i,:) = aver_xx;

end

miu_possibility = [];
estS.miu_aver = [];         %初始化一个空的IMM-ber滤波值
[p,q]=size(est_MC_miu_SIMM{MC,:});
for i=1:p        %3种模型
    xx = [];     %
    for  l=1:MC  %MC次循环
    xx1 = est_MC_miu_SIMM {l};
    xx = [xx ; xx1(i,:)];
    end
    aver_xx=sum(xx,1)/MC;         %貌似是求MC次滤波值的均值
    miu_possibility(i,:) = aver_xx;  %可能性

end
% 计算每列的和
column_sums = sum(miu_possibility, 1);
% 对每个元素进行归一化（每列元素除以该列的和）%可能性-->概率
estS.miu_aver = miu_possibility ./ column_sums; % 注意这里使用点除(./)进行逐元素操作

%%  %绘制OSPA误差
figure; 
subplot(3,1,1);
plot(1:meas.K,ospa_1,'r-', 'LineWidth', 1);
hold  on;
plot(1:meas.K,ospa_1SIMM,'b-', 'LineWidth', 1);
grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); 
legend('IMM-PHD','SIMM-PHD');
xlabel('Time');ylabel('OSPA Dist');

subplot(3,1,2);
plot(1:meas.K,ospa_2,'r-', 'LineWidth', 1);
hold on;
plot(1:meas.K,ospa_2SIMM,'b-', 'LineWidth', 1); 
grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); 
legend('IMM-PHD','SIMM-PHD');
xlabel('Time');ylabel('OSPA Loc');

subplot(3,1,3);
plot(1:meas.K,ospa_3,'r-', 'LineWidth', 1);
hold on;
plot(1:meas.K,ospa_3SIMM,'b-', 'LineWidth', 1); 
grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]);
legend('IMM-PHD','SIMM-PHD');
xlabel('Time');
ylabel('OSPA Card');

%% 绘制势估计plot cardinality
figure; cardinality= gcf; 
subplot(2,1,1); box on; hold on;
stairs(1:meas.K,truth.N,'k'); 
plot(1:meas.K,est.N_aver,'k.');
grid on;
legend(gca,'True','IMM Estimated');
set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 max(truth.N)+1]);
xlabel('Time'); ylabel('Cardinality');

subplot(2,1,2); box on; hold on;
stairs(1:meas.K,truth.N,'k'); 
plot(1:meas.K,estS.N_aver,'b.');
grid on;
legend(gca,'True','SIMM Estimated');
set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 max(truth.N)+1]);
xlabel('Time'); ylabel('Cardinality');

%% plot miu
%IMM
figure
t = 1:meas.K;
l1 = plot(t, est.miu_aver(1, :), 'r', 'LineWidth', 1.5);%CV
hold on;
l2 = plot(t, est.miu_aver(2, :), 'b', 'LineWidth', 1.5);%CT
hold on;
l3 = plot(t, est.miu_aver(3, :), 'k', 'LineWidth', 1.5);%CT2
hold on;
grid on;
% 获取坐标轴范围
XlabRang = get(gca,'xlim');  % 获取横坐标范围
XlabMmin = XlabRang(1);
XlabMmax = XlabRang(2);
% YlabRang = get(gca,'ylim');
YlabRang = [0,1];
YlabMmin = YlabRang(1);
YlabMmax = YlabRang(2);

% 绘制背景区域
FaceAlpha = 0.05;  % 背景
RecAcc = fill([0 0 model.tbirth(2) model.tbirth(2)],[YlabMmin YlabMmax YlabMmax YlabMmin],'r','FaceColor','b','EdgeColor','none','FaceAlpha',FaceAlpha);
RecC = fill([model.tbirth(2) model.tbirth(2) model.tdeath(2) model.tdeath(2)],[YlabMmin YlabMmax YlabMmax YlabMmin],'b','FaceColor','r','EdgeColor','none','FaceAlpha',FaceAlpha);
RecS = fill([model.tbirth(3)-1 model.tbirth(3)-1 model.tdeath(3) model.tdeath(3)],[YlabMmin YlabMmax YlabMmax YlabMmin],'r','FaceColor','k','EdgeColor','none','FaceAlpha',FaceAlpha);

% 最关键的 设置图像前后关系
set(gca,'child',[l1 l2 l3 RecAcc RecC RecS])  % 将曲线放在最前端

title('IMM-GMPHD Probabilities');
xlabel('Time (Moments)');
ylabel('Probability');
legend([l1,l2,l3],'CV Model', 'CT Model', 'CT Model');
hold off;

%SIMM
figure
t = 1:meas.K;
l4 = plot(t, estS.miu_aver(1, :), 'r', 'LineWidth', 1.5);
hold on;
l5 = plot(t, estS.miu_aver(2, :), 'b', 'LineWidth', 1.5);
hold on;
l6 = plot(t, estS.miu_aver(3, :), 'k', 'LineWidth', 1.5);
hold on;
grid on;
% 获取坐标轴范围
XlabRang = get(gca,'xlim');  % 获取横坐标范围
XlabMmin = XlabRang(1);
XlabMmax = XlabRang(2);
% YlabRang = get(gca,'ylim');
YlabRang = [0,1];
YlabMmin = YlabRang(1);
YlabMmax = YlabRang(2);

% 绘制背景区域
FaceAlpha = 0.05;  % 背景
RecAcc = fill([0 0 model.tbirth(2) model.tbirth(2)],[YlabMmin YlabMmax YlabMmax YlabMmin],'r','FaceColor','b','EdgeColor','none','FaceAlpha',FaceAlpha);
RecC = fill([model.tbirth(2) model.tbirth(2) model.tdeath(2) model.tdeath(2)],[YlabMmin YlabMmax YlabMmax YlabMmin],'r','FaceColor','r','EdgeColor','none','FaceAlpha',FaceAlpha);
RecS = fill([model.tbirth(3)-1 model.tbirth(3)-1 model.tdeath(3) model.tdeath(3)],[YlabMmin YlabMmax YlabMmax YlabMmin],'r','FaceColor','k','EdgeColor','none','FaceAlpha',FaceAlpha);

% 最关键的 设置图像前后关系
set(gca,'child',[l4 l5 l6 RecAcc RecC RecS])  % 将曲线放在最前端

title('SIMM-GMPHD Probabilities');
xlabel('Time (Moments)');
ylabel('Probability');
legend([l4,l5,l6],'CV Model', 'CT Model', 'CT Model');
hold off;


%% ————————————————————————————————————————————————————————————————
%% 初始化
est.error_aver=cell(truth.K,1);         
[p,q]=size(est_MC_error);
for i=1:q
    xx = [];
    for  l=1:MC
    xx1 = est_MC_error (l,i);
    xx = [xx , cell2mat(xx1)];
    end
    aver_xx=sum(xx,2)/MC;         %貌似是求MC次滤波值的均值
    est.error_aver{i} =[est.error_aver{i} aver_xx];

end

%%
estS.error_aver=cell(truth.K,1);         %初始化一个空的IMM-ber滤波值
[p_SIMM,q_SIMM]=size(est_MC_error_SIMM);
for i=1:q_SIMM
    xx = [];
    for  l=1:MC
    xx1 = est_MC_error_SIMM (l,i);
    xx = [xx , cell2mat(xx1)];
    end
    aver_xx_SIMM=sum(xx,2)/MC;         %貌似是求MC次滤波值的均值
    estS.error_aver{i} =[estS.error_aver{i} aver_xx_SIMM];

end

%%
%% 计算每个时间步的RMSE
%IMM
for k = 1:meas.K   
    error_matrix = est.error_aver{k}; % 提取当前时间步的误差矩阵
    squared_error = error_matrix.^2;% 计算每个状态变量的平方误差
    % 计算RMSE
    for i = 1:4
        IMM_rmse(k, i) = sqrt(mean(squared_error(i, :)));
    end
end
%SIMM
for k = 1:meas.K   
    error_matrix = estS.error_aver{k}; % 提取当前时间步的误差矩阵
    squared_error = error_matrix.^2;% 计算每个状态变量的平方误差
    % 计算RMSE
    for i = 1:4
        SIMM_rmse(k, i) = sqrt(mean(squared_error(i, :)));
    end
end

%% 绘制RMSE图
figure;
% 绘制X状态变量的RMSE曲线
plot(1:meas.K, IMM_rmse(:,1), 'LineWidth', 2, 'Color', 'r'); % x
hold on;
plot(1:meas.K, SIMM_rmse(:,1), 'LineWidth', 2, 'Color', 'b'); % x

% 添加图例和标签
% legend('x方向位置误差', 'x方向速度误差', 'y方向位置误差', 'y方向速度误差', 'Location', 'best');
legend('Estimation Error - IMM',  'Estimation Error - SIMM', 'Location', 'best');
xlabel('scans');
ylabel('meter');
title('Position RMSE in x-Axis');
grid on;
% 设置坐标轴范围
xlim([1, meas.K]);
ylim([0, max(max(IMM_rmse))*1.1]); % 稍微扩展y轴上限

figure;
% 绘制Y状态变量的RMSE曲线
plot(1:meas.K, IMM_rmse(:,3), 'LineWidth', 2, 'Color', 'r'); % y
hold on;
plot(1:meas.K, SIMM_rmse(:,3), 'LineWidth', 2, 'Color', 'b'); % y
% 添加图例和标签
% legend('x方向位置误差', 'x方向速度误差', 'y方向位置误差', 'y方向速度误差', 'Location', 'best');
legend('Estimation Error - IMM',  'Estimation Error - SIMM', 'Location', 'best');
xlabel('scans');
ylabel('meter');
title('Position RMSE in y-Axis');
grid on;
% 设置坐标轴范围
xlim([1, meas.K]);
ylim([0, max(max(IMM_rmse))*1.1]); % 稍微扩展y轴上限

% plot(1:K, rmse(:,2), 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250]); % vx
% plot(1:K, rmse(:,4), 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880]); % vy

toc

%子函数
function [X_track,k_birth,k_death]= extract_tracks(X,track_list,total_tracks)

K= size(X,1); 
x_dim= size(X{K},1); 
k=K-1; while x_dim==0, x_dim= size(X{k},1); k= k-1; end
X_track= zeros(x_dim,K,total_tracks);
k_birth= zeros(total_tracks,1);
k_death= zeros(total_tracks,1);

max_idx= 0;
for k=1:K
    if ~isempty(X{k})
        X_track(:,k,track_list{k})= X{k};
    end
    if max(track_list{k})> max_idx %new target born?
        idx= find(track_list{k}> max_idx);
        k_birth(track_list{k}(idx))= k;
    end
    if ~isempty(track_list{k}), max_idx= max(track_list{k}); 
    k_death(track_list{k})= k;
    end
end
end

function Xc= get_comps(X,c)
if isempty(X)
    Xc= [];       % 如果 X 为空矩阵，那么输出的 Xc 也为空矩阵
else
    Xc= X(c,:);   % 如果 X 不为空，根据索引向量 c 提取 X 中的行，将结果赋值给 Xc
                  % 这里的冒号 : 表示选取所有列，即提取出 X 中由 c 指定的行的所有列元素
end
end