function A=steering_matrix(N,theta)
        ob=(0:(N-1)).';
        A=exp(1j*ob*theta)/sqrt(N);
end