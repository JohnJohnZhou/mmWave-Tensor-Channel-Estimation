%{
Refrence:

''Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen and R. S. Blum, "Low-Rank 
Tensor Decomposition-Aided Channel Estimation for Millimeter Wave MIMO-OFDM 
Systems," in IEEE Journal on Selected Areas in Communications, 
vol. 35, no. 7, pp. 1524-1538, July 2017.''


Function: frequency extraction on steering vector
Date: Oct./2016
Author: Zhou Zhou

%}   

function f=extraction_f(N,R,Ob,dictionary_granularity)
    
    
            nn=N*dictionary_granularity;
            grid=((1:nn)*2*pi)/nn;
            D=exp(1j*(0:(N-1))'*grid);
            D=R*D;
            [~,index]=max(abs(normalization_column(D)'*normalization_column(Ob)));
            f=grid(index);
            
            
                function [A,lambda]=normalization_column(A)
                    lambda=sqrt(sum(abs(A).^2,1));
                    M=size(A,1);
                    A=A./(ones(M,1)*lambda);
                end
    end

