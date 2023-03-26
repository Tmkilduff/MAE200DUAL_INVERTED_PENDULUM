function [X] = RK4_step4(X00,E,A,B,K)

X=zeros(6,1000);
X(:,1) = X00;
h = 0.01;

for i = 1:1000

    f1=SSM(X(:,i),E,A,B,K); 
    f2=SSM(X(:,i)+f1/2*h,E,A,B,K); 
    f3=SSM(X(:,i)+f2/2*h,E,A,B,K); 
    f4=SSM(X(:,i)+h*f3,E,A,B,K); 

    X(:,i+1) = X(:,i) + h*(f1/6+(f2+f3)/3+f4/6); 
end

end