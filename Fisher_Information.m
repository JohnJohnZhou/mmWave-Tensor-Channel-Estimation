%{
Refrence:

''Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen and R. S. Blum, "Low-Rank 
Tensor Decomposition-Aided Channel Estimation for Millimeter Wave MIMO-OFDM 
Systems," in IEEE Journal on Selected Areas in Communications, 
vol. 35, no. 7, pp. 1524-1538, July 2017.''


Function: Fisher Information Matrix
Date: Oct./2016
Author: Zhou Zhou

%}

function [ Omega ] = Fisher_Information( theta_1,theta_2,theta_3,alpha,Q,P,S,sigma_2 )

L=4;

[M1,N1]=size(Q);
[M2,N2]=size(P);
[M3,N3]=size(S);


A=Q*steering_matrix(N1,theta_1);
B=P*steering_matrix(N2,theta_2);
C=S*steering_matrix(N3,theta_3)*diag(alpha);

A_tilde=Q*1j*diag(0:(N1-1))*steering_matrix(N1,theta_1);
B_tilde=P*1j*diag(0:(N2-1))*steering_matrix(N2,theta_2);
C_tilde=S*1j*diag(0:(N3-1))*steering_matrix(N3,theta_3)*diag(alpha);
G=S*steering_matrix(N3,theta_3);

num_matrix=1;
Omega_cell=cell.empty;
for mode_1=1:4
    for mode_2=1:4
        Omega_cell{num_matrix}=Covariance_Matrix( A,B,C,A_tilde,B_tilde,C_tilde,G,mode_1,mode_2,sigma_2);
        Omega_cell{num_matrix}=convert_covari_vector2matrix(Omega_cell{num_matrix},L);
        num_matrix=num_matrix+1;
    end
end

Omega_1=2*real([Omega_cell{1},Omega_cell{2},Omega_cell{3};...
    Omega_cell{5},Omega_cell{6},Omega_cell{7};...
    Omega_cell{9},Omega_cell{10},Omega_cell{11}]);
Omega_2=([Omega_cell{4};Omega_cell{8};Omega_cell{12}]);
Omega_3=([Omega_cell{13},Omega_cell{14},Omega_cell{15}]);
Omega_4=Omega_cell{16};
Omega=[Omega_1,Omega_2;Omega_3,Omega_4];

    
    
   

    function B=convert_covari_vector2matrix(A,L)
        % Convert the L^2*L^2 matrix A into L*L matrix B
        B=zeros(L,L);
        for l1=1:L
            for l2=1:L
                B(l1,l2)=A(L*(l1-1)+l1,L*(l2-1)+l2);
            end
        end
        
    end
end

