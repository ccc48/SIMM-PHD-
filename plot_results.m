function handles= plot_results(model,truth,meas,est,estSIMM)

[X_track,k_birth,k_death]= extract_tracks(truth.X,truth.track_list,truth.total_tracks);


%% plot ground truths
limit= [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) ];
figure; truths= gcf; hold on;
for i=1:truth.total_tracks
    Pt= X_track(:,k_birth(i):1:k_death(i),i); Pt=Pt([1 3],:);%x,y
    plot( Pt(1,:),Pt(2,:),'k-'); 
    plot( Pt(1,1), Pt(2,1), 'ko','MarkerSize',6);
    plot( Pt(1,31), Pt(2,31), 'ro','MarkerSize',6);
    plot( Pt(1,61), Pt(2,61), 'bo','MarkerSize',6);
    plot( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)), 'k^','MarkerSize',6);
end
axis equal; axis(limit); title('Ground Truths');
xlabel('x-coordinate (m)'); ylabel('y-coordinate (m)');


%%  plot tracks and measurements in x/y
figure; tracking= gcf; hold on;
subplot(211); box on; 

%plot x measurement
for k=1:meas.K
    if ~isempty(meas.Z{k})
        hlined= line(k*ones(size(meas.Z{k},2),1),meas.Z{k}(2,:).*cos(meas.Z{k}(1,:)),'LineStyle','none','Marker','x','Markersize',5,'Color',0.7*ones(1,3));
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
        hline2= line(k*ones(size(est.X_aver{k},2),1),P(1,:),'LineStyle','none','Marker','.','Markersize',8,'Color',0*ones(1,3));
    end
end

% plot x estimate(SIMM
for k=1:meas.K
    if ~isempty(estSIMM.X_aver{k})
        P= estSIMM.X_aver{k}([1 3],:);
        hline3= line(k*ones(size(estSIMM.X_aver{k},2),1),P(1,:),'LineStyle','none','Marker','.','Markersize',8,'Color','b');
    end
end
% 
subplot(212); box on;

%plot y measurement  
for k=1:meas.K
    if ~isempty(meas.Z{k})
        yhlined= line(k*ones(size(meas.Z{k},2),1),meas.Z{k}(2,:).*cos(meas.Z{k}(1,:)),'LineStyle','none','Marker','x','Markersize',5,'Color',0.7*ones(1,3));
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
        yhline2= line(k*ones(size(est.X_aver{k},2),1),P(2,:),'LineStyle','none','Marker','.','Markersize',8,'Color',0*ones(1,3));
    end
end

%plot y estimate(SIMM
for k=1:meas.K
    if ~isempty(estSIMM.X_aver{k})
        P= estSIMM.X_aver{k}([1 3],:);
        yhline3= line(k*ones(size(estSIMM.X_aver{k},2),1),P(2,:),'LineStyle','none','Marker','.','Markersize',8,'Color','b');
    end
end

subplot(211); xlabel('Time'); ylabel('x-coordinate (m)');
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(1,:));
legend([hline3 hline2 hline1 hlined],'SIMM Estimates','IMM Estimates','True tracks','Measurements');

subplot(212); xlabel('Time'); ylabel('y-coordinate (m)');
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(2,:));
% legend([yhline2 yhline1 yhlined],'Estimates          ','True tracks','Measurements');


% %plot cardinality
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
plot(1:meas.K,estSIMM.N_aver,'b.');
grid on;
legend(gca,'True','SIMM Estimated');
set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 max(truth.N)+1]);
xlabel('Time'); ylabel('Cardinality');
% return
% handles=[ truths tracking ospa cardinality ];
handles= [truths  tracking cardinality ];


function [X_track,k_birth,k_death]= extract_tracks(X,track_list,total_tracks)

K= size(X,1); 
x_dim= size(X{K},1); 
k=K-1; while x_dim==0, x_dim= size(X{k},1); k= k-1; end;
X_track= zeros(x_dim,K,total_tracks);
k_birth= zeros(total_tracks,1);
k_death= zeros(total_tracks,1);

max_idx= 0;
for k=1:K
    if ~isempty(X{k}),
        X_track(:,k,track_list{k})= X{k};
    end;
    if max(track_list{k})> max_idx, %new target born?
        idx= find(track_list{k}> max_idx);
        k_birth(track_list{k}(idx))= k;
    end;
    if ~isempty(track_list{k}), max_idx= max(track_list{k}); end;
    k_death(track_list{k})= k;
end;

function Xc= get_comps(X,c)

if isempty(X)
    Xc= [];
else
    Xc= X(c,:);
end
