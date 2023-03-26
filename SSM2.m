function dx_hat = SSM2(x_hat,A,B,C,L,K,X)
dx_hat = (A+L*C)*x_hat + (B*K-L*C)*X;
end