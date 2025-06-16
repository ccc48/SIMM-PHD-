function [H,U]= CT_ekf_update_mat(model,mu)

 p= mu([1 3],:);
 mag= p(1)^2 + p(2)^2;
 sqrt_mag= sqrt(mag);
 H= [ p(2)/mag  0  -p(1)/mag  0  0; ...
         p(1)/sqrt_mag  0  p(2)/sqrt_mag 0  0 ];

 % fprintf('p矩阵:\n');
 % disp(p);
 % fprintf('mag矩阵:\n');
 % disp(mag);
 % fprintf('sqrt_mag矩阵:\n');
 % disp(sqrt_mag);
 % fprintf('H矩阵:\n');
 % disp(H)
 
 U= eye(2);
end