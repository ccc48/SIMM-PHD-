function Z= CT_gen_observation_fn(model,X,W)

%r/t observation equation

if ~isnumeric(W)
    if strcmp(W,'noise')
        W= model.D_CT*randn(size(model.D_CT,2),size(X,2));
    elseif strcmp(W,'noiseless')
        W= zeros(size(model.D_CT,1),size(X,2));
    end
end

if isempty(X)
    Z= [];
else %modify below here for user specified measurement model
    % P= X([1 3],:);
    % Z(1,:)= atan2(P(1,:),P(2,:));   
    % Z(2,:)= sqrt(sum(P.^2));
    % Z= Z+ W;

    Z_1i = model.H_CT*X;  %ONE METHOD
    % delete omega
    Z_2i = [Z_1i([1 2],:)];
    %ADD NOISE
    Z = Z_2i + W;
end