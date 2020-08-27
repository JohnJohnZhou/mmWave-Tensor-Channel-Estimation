%{

Reference:

''Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen and R. S. Blum, "Low-Rank 
Tensor Decomposition-Aided Channel Estimation for Millimeter Wave MIMO-OFDM 
Systems," in IEEE Journal on Selected Areas in Communications, 
vol. 35, no. 7, pp. 1524-1538, July 2017.''

Function: the CRB and algorithm versus number of frames(M2)
Date: Oct./2016
Author: Zhou Zhou

%}

clc;
clear;
addpath('./tensor_toolbox_2.6')

M1=6;
N1=32;
M2=[2,3,4:4:64];
N2=64;
M3=2;
N3=128;
SNR=10;

num=1;


load('.\CSI_file\alpha.mat');
load('.\CSI_file\theta_1.mat');
load('.\CSI_file\theta_2.mat');
load('.\CSI_file\theta_3.mat');
load('.\CSI_file\L.mat')

crb_theta_1=zeros(length(M3),num);
crb_theta_2=zeros(length(M3),num);
crb_theta_3=zeros(length(M3),num);
crb_alpha=zeros(length(M3),num);
error_theta_1=zeros(length(M3),num);
error_theta_2=zeros(length(M3),num);
error_theta_3=zeros(length(M3),num);
error_alpha=zeros(length(M3),num);
error_H=zeros(length(M3),num);
ff=1;
Q=exp(1j*2*pi*rand(M1,N1));
P_all=exp(1j*2*pi*rand(M2(end),N2));
I=randperm(N3);
I=I(1:M3);
S=eye(N3);
S=S(I,:);
for mm=M2
mm
        P=P_all(1:mm,:);
    for nn=1:num

    %% Parameters
%         L=4;
%         theta_1=rand(1,L)*2*pi;
%         theta_2=rand(1,L)*2*pi;
%         theta_3=rand(1,L)*2*pi;
%         alpha=(randn(L,1)+1j*randn(L,1))/sqrt(2);

        [Y,H,sigma_2]=Observation(theta_1,theta_2,theta_3,alpha,Q,P,S,SNR);


        %% CRB
        CRB_value=CRB(theta_1,theta_2,theta_3,alpha,Q,P,S,sigma_2);
        CRB_value=diag(CRB_value);

        crb_theta_1(ff,nn)=sum(CRB_value(1:L));
        crb_theta_2(ff,nn)=sum(CRB_value((L+1):2*L));
        crb_theta_3(ff,nn)=sum(CRB_value((2*L+1):3*L));
        crb_alpha(ff,nn)=sum(CRB_value((3*L+1):4*L));

        %% CP based estimation
        dim_ratio=300;

        [est_H,est_alpha,est_theta_1,est_theta_2,est_theta_3,iter(ff,nn)]=Tensor_Parameters_Estimation(Y,Q,P,S,L,dim_ratio,0);

        [error_theta_1_,error_theta_2_,error_theta_3_,error_alpha_,error_H_]=parameter_result(est_H,est_alpha,est_theta_1,est_theta_2,est_theta_3,H,alpha,theta_1,theta_2,theta_3);

        error_theta_1(ff,nn)=error_theta_1_;

        error_theta_2(ff,nn)=error_theta_2_;

        error_theta_3(ff,nn)=error_theta_3_;

        error_alpha(ff,nn)=error_alpha_;

        error_H(ff,nn)=error_H_;
    end
     ff=ff+1;
end
tensor_Frame_MSE.M2=M2;
tensor_Frame_MSE.crb_theta_1=crb_theta_1;
tensor_Frame_MSE.crb_theta_2=crb_theta_2;
tensor_Frame_MSE.crb_theta_3=crb_theta_3;
tensor_Frame_MSE.crb_alpha=crb_alpha;
tensor_Frame_MSE.error_theta_1=error_theta_1;
tensor_Frame_MSE.error_theta_2=error_theta_2;
tensor_Frame_MSE.error_theta_3=error_theta_3;
tensor_Frame_MSE.error_alpha=error_alpha;
tensor_Frame_MSE.error_H=error_H;
tensor_Frame_MSE.iter=iter;




