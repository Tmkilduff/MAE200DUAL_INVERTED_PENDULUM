function dX = DRE_X(X,E,A,B,R,Q)
% dX = -E'^(-1)*A'*X - X*A*E^(-1) + X*B*R^(-1)*B'*X - E'^(-1)*Q*E^(-1);
dX = A'*X*E + E'*X*A - E'*X*B*R^(-1)*B'*X*E + Q;
dX = -E'^(-1)*dX*E^(-1);
end