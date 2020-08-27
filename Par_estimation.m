%{
Refrence:

''Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen and R. S. Blum, "Low-Rank 
Tensor Decomposition-Aided Channel Estimation for Millimeter Wave MIMO-OFDM 
Systems," in IEEE Journal on Selected Areas in Communications, 
vol. 35, no. 7, pp. 1524-1538, July 2017.''


Function: parameters estimation by factor matrix
Date: Oct./2016
Author: Zhou Zhou

%}

function [ est_mat_MS,est_mat_BS,est_mat_F,est_lambda,est_f_MS,est_f_BS,est_f_F ] = Par_estimation(Z,A,B,C,W,P,S,dim_ratio,sigma_2 )
     
% Parameter extraction
     N_MS=size(W,2);
     N_BS=size(P,2);
     F=size(S,2);
     L=size(A,2);
     [est_f_MS] = extraction_f(N_MS,W,conj(A),dim_ratio);
     [est_f_BS] = extraction_f(N_BS,P,conj(B),dim_ratio);
     [est_f_F] = extraction_f(F,S,conj(C),150);
     
     est_mat_MS=steering_matrix(N_MS,est_f_MS);
     est_mat_BS=steering_matrix(N_BS,est_f_BS);
     est_mat_F=steering_matrix(F,est_f_F);
     
     % first method: good
     Equal_A=khatrirao(P*est_mat_BS,W*est_mat_MS);
     Equal_y=tenmat(Z,3);% we only estimate the channel without pathloss
     Equal_y=Equal_y.data;
     
     est_C=(Equal_y)*pinv(Equal_A.');
     %est_C=(Equal_y)*conj(Equal_A)/(Equal_A.'*conj(Equal_A)+sigma_2*eye(L));
     est_lambda=scale(est_C,S*est_mat_F);    
     
% second method
     
%      lambda_1=scale(A,W*est_mat_MS);
%      lambda_2=scale(B,P*est_mat_BS);
%      lambda_3=scale(C,S*est_mat_F);
%      est_lambda=lambda_1.*lambda_2.*lambda_3;
 
     
    function f=extraction_f(N,R,Ob,dictionary_granularity)
            nn=N*dictionary_granularity;
            grid=((1:nn)*2*pi)/nn;
            D=exp(1j*(0:(N-1))'*grid);
            D=R*D;
            [~,index]=max(abs(normalization_column(D)'*normalization_column(Ob)));
            f=grid(index);
    end

    function [A,lambda]=normalization_column(A)
        lambda=sqrt(sum(abs(A).^2,1));
        M=size(A,1);
        A=A./(ones(M,1)*lambda);
    end

    function [H]=steering_matrix(N,theta)
        
        H=exp(1j*(0:(N-1))'*theta)/sqrt(N);

    end

    function lambda=scale(A,B)
        D=pinv(B)*A;
        lambda=diag(D);
    end
end

