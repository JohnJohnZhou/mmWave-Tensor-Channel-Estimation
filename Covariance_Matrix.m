%{
Refrence:

''Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen and R. S. Blum, "Low-Rank 
Tensor Decomposition-Aided Channel Estimation for Millimeter Wave MIMO-OFDM 
Systems," in IEEE Journal on Selected Areas in Communications, 
vol. 35, no. 7, pp. 1524-1538, July 2017.''


Function: Covariance Matrix Calculation fo Fisher Information Matrix
Date: Oct./2016
Author: Zhou Zhou

%}

function [ Cov ] = Covariance_Matrix( A,B,C,A_tilde,B_tilde,C_tilde,G,mode_1,mode_2,sigma_2 )

I=size(A,1);
J=size(B,1);
K=size(C,1);

if mode_1==1
    temp1=kron(A_tilde,khatrirao(C,B));
elseif mode_1==2
    temp1=kron(B_tilde,khatrirao(C,A));    
elseif mode_1==3
    temp1=kron(C_tilde,khatrirao(B,A));
elseif mode_1==4
    temp1=kron(G,khatrirao(B,A));    
end

if mode_2==1
    temp2=kron(A_tilde,khatrirao(C,B));
elseif mode_2==2
    temp2=kron(B_tilde,khatrirao(C,A));    
elseif mode_2==3
    temp2=kron(C_tilde,khatrirao(B,A));
elseif mode_2==4
    temp2=kron(G,khatrirao(B,A));    
end

C_N=Noise_Covariance_Matrix( I,J,K,mode_1,mode_2,sigma_2 );
Cov=(temp1.'*conj(C_N)*conj(temp2))/(sigma_2^2);


end

