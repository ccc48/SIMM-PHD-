function [w_new,x_new,P_new]= gaus_merge_newnew(w,x,P,threshold)

L= length(w); x_dim= size(x,1);
I= 1:L;      %从 1 到 L 的索引数组 I，用于在后续循环中迭代高斯分量
el= 1;      %计数器 el，用于跟踪合并后的新高斯分量的索引。

if all(w==0)
    w_new = [];
    x_new = [];
    P_new = [];
    return;
end
    w_new = [];
    x_new = [];
    P_new = [];

while ~isempty(I),
    [notused,j]= max(w); j= j(1);     %找到权重数组 w 中的最大值notused及其索引 j
    Ij= []; iPt= inv(P(:,:,j));      %使用 inv 函数可能导致数值不稳定，从而出现警告：矩阵在工作精度内为奇异的

    for i= I   %遍历剩余的高斯分量索引 
        val= (x(:,i)-x(:,j))'*iPt*(x(:,i)-x(:,j));%计算每个高斯分量与当前最大权重高斯分量 j 之间的马氏距离
        if val <= threshold,  %如果距离小于或等于阈值 threshold
            Ij= [ Ij i ];    %则将该高斯分量的索引添加到 Ij
        end;
    end;

    %判断Ij是否非空
    if isempty(Ij)
        w_new = [w_new; w(j)];    %若Ij为空矩阵，则保持wmP不变，即没有可以合并的高斯分量
        x_new = [x_new x(:,j)];
        P_new = cat(3,P_new, P(:,:,j));

        % x_new(:,el)= x_new(:,el)/w_new(el);      %没有合并，不需要对均值和协方差矩阵进行标准化
        % P_new(:,:,el)= P_new(:,:,el)/w_new(el);
        I= setdiff(I,j);       %这行更新了索引数组 I，移除了已经合并的高斯分量的索引。
        w(j)= -1;
        el= el+1;

     
    else              %若Ij非空，则计算合并后的新高斯分量
   w_new= [w_new; sum(w(Ij))];    %计算了合并后的新高斯分量de
   x_new= [x_new wsumvec(w(Ij),x(:,Ij),x_dim)];
   P_new= cat(3,P_new, wsummat(w(Ij),P(:,:,Ij),x_dim));

   x_new(:,el)= x_new(:,el)/w_new(el);      %对合并后的均值和协方差矩阵进行了标准化
   P_new(:,:,el)= P_new(:,:,el)/w_new(el);
   I= setdiff(I,Ij);       %这行更新了索引数组 I，移除了已经合并的高斯分量的索引。
   w(Ij)= -1;
   el= el+1;
    end
end

function out = wsumvec(w,vecstack,xdim)
    wmat = repmat(w',[xdim,1]);
    out  = sum(wmat.*vecstack,2);

function out = wsummat(w,matstack,xdim)
    w = reshape(w,[1,1,size(w)]);
    wmat = repmat(w,[xdim,xdim,1]);
    out = sum(wmat.*matstack,3);
    
