%{
Refrence:

''Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen and R. S. Blum, "Low-Rank 
Tensor Decomposition-Aided Channel Estimation for Millimeter Wave MIMO-OFDM 
Systems," in IEEE Journal on Selected Areas in Communications, 
vol. 35, no. 7, pp. 1524-1538, July 2017.''


Function: ALS for CP decomposition
Date: Oct./2016
Author: Zhou Zhou

%}

function [ X,A,B,C,res,initer ] = Complex_ALS(Z,rank,sigma_2 )
% 
    I=size(Z);
    maxiter=3e3;

    A=tenmat(Z,1);
    A=take_svd(A.data,rank);
    
    B=tenmat(Z,2);
    B=take_svd(B.data,rank);

    C=tenmat(Z,3);
    C=take_svd(C.data,rank);
%     rng(20)
%     A=randn(I(1),rank)+1j*randn(I(1),rank);
%     rng(20)
%     B=randn(I(2),rank)+1j*randn(I(2),rank);
%     rng(20)
%     C=randn(I(3),rank)+1j*randn(I(3),rank);

    terminated=0;
    outiter=0;
    initer=0;
    
    X=ktensor(ones(rank,1),A,B,C);
    X=tensor(X);
    
    global verbose
    
    while ~terminated
        
        old_X=X;
            
        
        A=updatafactor(1,khatrirao(C,B));
        [A,~]=normalization_column(A);
        B=updatafactor(2,khatrirao(C,A)); 
        [B,~]=normalization_column(B);
        C=updatafactor(3,khatrirao(B,A));
        [C,lambda]=normalization_column(C);
        X=ktensor(lambda.',A,B,C);
        X=tensor(X);

        initer=initer+1;
        
        res=norm(Z-X)/norm(Z);
        
        
        if norm(X-old_X)/norm(X)<1e-8||outiter>maxiter
            if outiter<=maxiter
                converge=1;
            else
                converge=0;
            end
            terminated=1;
            C=C*diag(lambda);
        end
        
        outiter=outiter+1;
        
        if verbose==1
            fprintf('The number of iteration is %6.2f\n',outiter);
            fprintf('The residuals are %2.10f\n',norm(Z-X)/norm(Z));
            if  terminated==1
                  fprintf('\n\n');
            end
        end
        
    end

    
    
    function A=updatafactor(n_mode,PI)
        Z_fold=tenmat(Z,n_mode);
        nn=size(PI,2);
        new_A=(conj(double(Z_fold))*PI)/(PI'*PI+sigma_2*eye(nn));
        A=conj(new_A);
    end

    function U_matr=take_svd(Matr,rank)
        [U_matr,S,~]=svd(Matr,'econ');
        [m,n]=size(U_matr);
        if n>=rank
            U_matr=U_matr(:,1:rank);%*sqrt(S(1:rank,1:rank))
        else
            rng(10,'v5uniform');
            U_matr=[U_matr,randn(m,rank-n)+1j*randn(m,rank-n)];
        end
    end

    function [A,lambda]=normalization_column(A)
        lambda=sqrt(sum(abs(A).^2,1));
        M=size(A,1);
        A=A./(ones(M,1)*lambda);
    end

end

