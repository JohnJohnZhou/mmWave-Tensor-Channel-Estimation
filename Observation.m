function [ Y, H, sigma_2 ] = Observation( theta_1,theta_2,theta_3,alpha,Q,P,S,SNR )

[M1,N1]=size(Q);
[M2,N2]=size(P);
[M3,N3]=size(S);

A=steering_matrix(N1,theta_1);
B=steering_matrix(N2,theta_2);
C=steering_matrix(N3,theta_3);
Y=ktensor(alpha,Q*A,P*B,S*C);
Y=tensor(Y);

H=ktensor(alpha,A,B,C);
H=tensor(H);

N=randn(M1,M2,M3)+1j*randn(M1,M2,M3);
N=N/sqrt(2);
N=tensor(N);

sigma=(norm(Y)/norm(N))*10^(-SNR/20);
sigma_2=(sigma)^2;

N=N*sigma;

Y=Y+N;


end

