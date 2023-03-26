function [P] = RK4_forwards(P0,E,A,C,Q1,Q2)
P=zeros(size(A,1),size(A,2),size(A,3));
P(:,:,end) = P0;
h = 0.01;

for i = 1:size(P,3):2 %size(X,3):-1:2
    f1=DRE_P(P(:,:,i), E(:,:,i), A(:,:,i), C, Q1, Q2); 
    f2=DRE_P(P(:,:,i)+f1/2*h, E(:,:,i), A(:,:,i), C, Q1, Q2); 
    f3=DRE_P(P(:,:,i)+f2/2*h, E(:,:,i), A(:,:,i), C, Q1, Q2); 
    f4=DRE_P(P(:,:,i)+h*f3, E(:,:,i), A(:,:,i), C, Q1, Q2); 

    P(:,:,i+1) = P(:,:,i) + h*(f1/6+(f2+f3)/3+f4/6); 
end
end