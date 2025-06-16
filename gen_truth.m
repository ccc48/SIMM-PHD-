function truth= gen_truth(model)
%variables
truth.K= 150;                   %length of data/number of scans
truth.X= cell(truth.K,1);             %ground truth for states of targets  
truth.N= zeros(truth.K,1);            %ground truth for number of targets
truth.L= cell(truth.K,1);             %ground truth for labels of targets (k,i)
truth.track_list= cell(truth.K,1);    %absolute index target identities (plotting)
truth.total_tracks= 0;          %total number of appearing tracks

%target initial states and birth/death times
nbirths= 1;     %目标数
                %x  vx  y  vy 
% xstart(:,1)  = [ 1300; 25;  1200; 0];       tbirth(1)  = 1;     tdeath(1)  = truth.K+1;
% xstart(:,1)  = [ 13000; 250;  12000; 0];       tbirth(1)  = 1;     tdeath(1)  = truth.K+1;
xstart(:,1)  = [ 2100; 40;  2000; 0];       tbirth(1)  = 1;     tdeath(1)  = truth.K+1;
% xstart(:,1)  = [ 2100; 40;  2800; 0];       tbirth(1)  = 1;     tdeath(1)  = truth.K+1;

% start with CT1
G = model.G{2};  % Start with CT1
Q = model.Q{2};

for targetnum=1:nbirths
    for k = 1:truth.K
        if any(model.maneuvers(1,:) == k)
            model.idx = model.maneuvers(2, model.maneuvers(1,:) == k);  %maneuvers(1,:) == k：这是一个逻辑表达式，用于比较maneuvers矩阵第一行的每个元素是否等于k。结果是一个逻辑数组，其中等于k的位置为true，其他位置为false。即[ 1  0  0]
            switch model.idx
                case 1   % % 如果model_idx为1  ,即CT
                    G = model.G{2};  % Start with CT
                    Q = model.Q{2};
                    truth.modelType = 'CT';
                    targetstate = xstart;
                    for k=model.tbirth(model.idx):min(model.tdeath(model.idx),truth.K)
                        targetstate = gen_newstate_fn(k,model,targetstate,G,Q,'noiseless');
                        truth.X{k}= [truth.X{k} targetstate];
                        truth.track_list{k} = [truth.track_list{k} targetnum];
                        truth.N(k) = truth.N(k) + 1;
                    end
                case 2   %CV
                    G = model.G{1};  
                    Q = model.Q{1};
                    truth.modelType = 'CV';
                    targetstate = targetstate;
                    for k=model.tbirth(model.idx):min(model.tdeath(model.idx),truth.K)
                        targetstate = gen_newstate_fn(k,model,targetstate,G,Q,'noiseless');
                        truth.X{k}= [truth.X{k} targetstate];
                        truth.track_list{k} = [truth.track_list{k} targetnum];
                        truth.N(k) = truth.N(k) + 1;
                    end
                    
                case 3   %CT
                    G = model.G{3};  
                    Q = model.Q{3};
                    truth.modelType = 'CT';
                    targetstate = targetstate;
                    for k=model.tbirth(model.idx):min(model.tdeath(model.idx),truth.K)
                        targetstate = gen_newstate_fn(k,model,targetstate,G,Q,'noiseless');
                        truth.X{k}= [truth.X{k} targetstate];
                        truth.track_list{k} = [truth.track_list{k} targetnum];
                        truth.N(k) = truth.N(k) + 1;
                    end
            end
        end
  
    truth.total_tracks= nbirths;
    end
    end
