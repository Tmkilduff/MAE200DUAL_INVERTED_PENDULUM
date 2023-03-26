function [X] = RK4_backwards(X0,E,A,B,R,Q)
X=zeros(size(A,1),size(A,2),size(A,3));
X(:,:,end) = X0;
h = 0.01;

for i = size(X,3):-1:2
    f1=DRE_X(X(:,:,i), E(:,:,i), A(:,:,i), B, R, Q); 
    f2=DRE_X(X(:,:,i)-f1/2*h, E(:,:,i), A(:,:,i), B, R, Q); 
    f3=DRE_X(X(:,:,i)-f2/2*h, E(:,:,i), A(:,:,i), B, R, Q); 
    f4=DRE_X(X(:,:,i)-h*f3, E(:,:,i), A(:,:,i), B, R, Q); 

    X(:,:,i-1) = X(:,:,i) - h*(f1/6+(f2+f3)/3+f4/6); 
end

end