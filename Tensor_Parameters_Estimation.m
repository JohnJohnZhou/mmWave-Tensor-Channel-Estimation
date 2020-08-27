%{
Refrence:

''Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen and R. S. Blum, "Low-Rank 
Tensor Decomposition-Aided Channel Estimation for Millimeter Wave MIMO-OFDM 
Systems," in IEEE Journal on Selected Areas in Communications, 
vol. 35, no. 7, pp. 1524-1538, July 2017.''


Function: Tensor decomposition aided method for parameters estimation
Date: Oct./2016
Author: Zhou Zhou

%}

function [est_H,alpha,theta_MS,theta_BS,tao,initer]=Tensor_Parameters_Estimation(Y,W,P,S,L,dim_ratio,sigma_2)


[est_Y,F_MS,F_BS,F_S,res,initer]=Complex_ALS(Y,L,sigma_2);

[est_f_MS,est_f_BS,est_f_F,alpha,theta_MS,theta_BS,tao]=Par_estimation(est_Y,F_MS,F_BS,F_S,W,P,S,dim_ratio,sigma_2);

est_H=ktensor(alpha,est_f_MS,est_f_BS,est_f_F);

est_H=tensor(est_H);

end