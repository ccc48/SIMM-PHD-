function [w_update,m_update,P_update,L_update,LK] = CVfilter(k,model,meas,est,filter,w_update,m_update,P_update,L_update)
%=== Filtering 
    %---prediction 
    [m_predict,P_predict] = CV_kalman_predict(model,m_update,P_update);                       %surviving components
    w_predict= model.P_S*w_update;                                                                  %surviving weights

    m_predict= cat(2,model.m_birth,m_predict); P_predict=cat(3,model.P_birth,P_predict);            %append birth components
    w_predict= cat(1,model.w_birth,w_predict)  ;                                                   %append birth weights
                                                 
    L_predict= model.L_birth+L_update;                                                              %number of predicted components
    
    %---gating
    if filter.gate_flag
        meas.Z{k} = gate_meas_gms(meas.Z{k},filter.gamma,model,m_predict,P_predict);        
    end
        
    %---update
    %number of measurements
    m= size(meas.Z{k},2);
    
    %missed detection term 
    w_update = model.Q_D*w_predict;
    m_update = m_predict;
    P_update = P_predict;
    
    if m~=0
        %m detection terms 
        [qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k},model,m_predict,P_predict);
        for ell=1:m
            w_temp = model.P_D*w_predict(:).*qz_temp(:,ell); 
            w_temp= w_temp./(model.lambda_c*model.pdf_c + sum(w_temp));
            w_update = cat(1,w_update,w_temp);
            m_update = cat(2,m_update,m_temp(:,:,ell));
            P_update = cat(3,P_update,P_temp);
        end
    end         
            
    %Likelihood
    Z = meas.Z{k};
    LK = CV_compute_likelihood(k,model,Z,est,m_predict,P_predict);
    


    %---mixture management
    L_posterior= length(w_update);
    
    %pruning, merging, capping
    [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,filter.elim_threshold);    L_prune= length(w_update);
    [w_update,m_update,P_update]= gaus_merge_newnew(w_update,m_update,P_update,filter.merge_threshold);   L_merge= length(w_update);
    [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,filter.L_max);
    L_cap  = length(w_update);
    
    L_update= L_cap;

    %--- state extraction
    idx= find(w_update > 0.5 );
    P1 = [];
    for j=1:length(idx)
        repeat_num_targets= round(w_update(idx(j)));
        est.W{k,1}= [ est.W{k,1}; repmat(w_update(idx(j),:),[repeat_num_targets,1]) ];
        est.X{k,1}= [ est.X{k,1} repmat(m_update(:,idx(j)),[1,repeat_num_targets]) ];
        for rll = 1:repeat_num_targets
            P = P_update(:,:,idx(j)) ;
            P1 = cat(3,P1,P);
        end

        est.N(k,1)= est.N(k,1)+repeat_num_targets;
        est.L{k,1}= [];
    end
    est.P{k,1} = P1;
end
    
 

function LK = CV_compute_likelihood(k,model,Z,est,m_predict,P_predict)

    n = model.z_dim;
    num_GM = size(m_predict,2);

    if isempty(m_predict)|isempty(Z)
        LK= [];
    else

    %残差  gating后，量测为空则残差为空
    numZ = size(Z,2);
    LK=[];
    for Gll = 1:num_GM
        v = Z -model.H*repmat(m_predict(:,Gll),[1 numZ]);
        %残差为空则似然为空
        %残差协方差
        S = model.H *P_predict(:,:,Gll)*model.H' + model.R;
        Vs= chol(S);        %Cholesky分解
        det_S= prod(diag(Vs))^2; %S的行列式
        inv_sqrt_S= inv(Vs);    %S的逆的平方根
        iS= inv_sqrt_S*inv_sqrt_S';    %S的逆
        %残差为空，则似然为空
        for vll = 1:size(v,2)
            Likel = exp(-0.5*n*log(2*pi)-0.5*log(det_S)-0.5*v(:,vll)'*iS*v(:,vll));  %似然
                
            %存
            LK= [LK  Likel];
            
        end
    end
    end 
end