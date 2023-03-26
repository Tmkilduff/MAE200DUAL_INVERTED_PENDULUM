function [x_k] = RungeKutta4(x,u_k,s)

for n=1:s.N
    u=u_k(n);       % (from t=0 -> T), compute cost
    f1=RHS(x,u,s); 
    f2=RHS(x+s.h*f1/2,u,s); 
    f3=RHS(x+s.h*f2/2,u,s); 
    f4=RHS(x+s.h*f3,u,s);

    x=x+s.h*(f1/6+(f2+f3)/3+f4/6); 
    x_k(1:6,n+1)=x; 
    u=u_k(n+1);
    x_k(7:9,n)=f1(4:6); 
end
%     if n==s.N 
%         c=.25; 
%     end
%     J=J+c*s.h*(x'*s.Q*x+u'*s.R*u);
end