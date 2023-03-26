function dX = SSM(X,E,A,B,K)
dX = E\(A+B*K)*X;
end