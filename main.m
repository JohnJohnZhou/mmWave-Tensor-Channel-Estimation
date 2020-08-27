%{
Refrence:

''Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen and R. S. Blum, "Low-Rank 
Tensor Decomposition-Aided Channel Estimation for Millimeter Wave MIMO-OFDM 
Systems," in IEEE Journal on Selected Areas in Communications, 
vol. 35, no. 7, pp. 1524-1538, July 2017.''


Function: Main testing
Date: Oct./2016
Author: Zhou Zhou

%}

clc;
clear;
addpath('./tensor_toolbox_2.6')
%% Parameters
load('.\CSI_file\alpha.mat');
load('.\CSI_file\theta_1.mat');
load('.\CSI_file\theta_2.mat');
load('.\CSI_file\theta_3.mat');
load('.\CSI_file\L.mat')

M1=6;
N1=32;
M2=6;
N2=64;
M3=6;
N3=128;
SNR=30;

Q=exp(1j*2*pi*rand(M1,N1))/sqrt(N1);
P=exp(1j*2*pi*rand(M2,N2))/sqrt(N2);
I=randperm(N3);
I=I(1:M3);
S=eye(N3);
S=S(I,:);

[Y,H,sigma_2]=Observation(theta_1,theta_2,theta_3,alpha,Q,P,S,SNR);


%% CRB
CRB_value=CRB(theta_1,theta_2,theta_3,alpha,Q,P,S,sigma_2);
CRB_value=diag(CRB_value);

crb_tehta1 = real(sum(CRB_value(1:L)))
crb_tehta2 = real(sum(CRB_value((L+1):2*L)))
crb_tehta3 = real(sum(CRB_value((2*L+1):3*L)))
crb_alpha = real(sum(CRB_value((3*L+1):4*L)))

%% CP based estimation
dim_ratio=200;

[est_H,est_alpha,est_theta_1,est_theta_2,est_theta_3]=Tensor_Parameters_Estimation(Y,Q,P,S,L,dim_ratio,0);

[error_theta_1,error_theta_2,error_theta_3,error_alpha,error_H]=parameter_result(est_H,est_alpha,est_theta_1,est_theta_2,est_theta_3,H,alpha,theta_1,theta_2,theta_3);

error_theta_1

error_theta_2

error_theta_3

error_alpha

error_H