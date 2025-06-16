function [F,G]= ekf_predict_mat_CA(model,mu_old)

tol= 1e-6;
T= model.T;

A = [ 1 T; 0 1];
AA = [(T^2)/2; T];
AAA = [1 ; 1];
F= [ A zeros(2,2) AA zeros(2,1);
     zeros(2,2)   A  zeros(2,1) AA;
     zeros(2,5)  AAA];

G= model.B_CA;
    