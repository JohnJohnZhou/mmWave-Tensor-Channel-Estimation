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

function [ C ] = Noise_Covariance_Matrix( I,J,K,mode_1,mode_2,sigma_2 )
% The cross expectation matrix of Circular Symmetric Gaussian
% Input: mode_1: 1--A;2--B;3--C
%        sigma_2: the variance of noise


if mode_1==mode_2
    C=eye(I*J*K);
else
    C=zeros(I*J*K);
    for i=1:I
        for j=1:J
            for k=1:K
                
                if mode_1==1
                    index_1=j+(k-1)*J+(i-1)*J*K;
                elseif mode_1==2
                    index_1=i+(k-1)*I+(j-1)*I*K;
                elseif mode_1==3
                    index_1=i+(j-1)*I+(k-1)*I*J;
                elseif mode_1==4
                    index_1=i+(j-1)*I+(k-1)*I*J;
                end
                
                if mode_2==1
                    index_2=j+(k-1)*J+(i-1)*J*K;
                elseif mode_2==2
                    index_2=i+(k-1)*I+(j-1)*I*K;
                elseif mode_2==3
                    index_2=i+(j-1)*I+(k-1)*I*J;
                elseif mode_2==4
                    index_2=i+(j-1)*I+(k-1)*I*J;                    
                end
                    
                C(index_1,index_2)=sigma_2;
                
            end
        end
    end
end





end

