function x = CartoPol(m)
    P= m([1 2],:);
    x(1,:)= atan2(P(1,:),P(2,:));   
    x(2,:)= sqrt(sum(P.^2));
end