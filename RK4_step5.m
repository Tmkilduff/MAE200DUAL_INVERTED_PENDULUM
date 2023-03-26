function [x_hat] = RK4_step5(X0,A,B,C,L,K,X_new)

x_hat=zeros(6,1000);
x_hat(:,1) = X0;
h = 0.01;

for i = 1:1000
    f1=SSM2(x_hat(:,i),A,B,C,L,K,X_new(:,i)); 
    f2=SSM2(x_hat(:,i)+f1/2*h, A,B,C,L,K,X_new(:,i)); 
    f3=SSM2(x_hat(:,i)+f2/2*h, A,B,C,L,K,X_new(:,i)); 
    f4=SSM2(x_hat(:,i)+h*f3, A,B,C,L,K,X_new(:,i)); 

    x_hat(:,i+1) = x_hat(:,i) + h*(f1/6+(f2+f3)/3+f4/6); 
end
end