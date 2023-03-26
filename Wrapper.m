%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tim Kilduff
%MAE 200
%Dual inverted pendulum
%https://github.com/Tmkilduff/MAE200DUAL_INVERTED_PENDULUM.git

clear all, close all, clc
% constants
T = 3; %sec
s.h=0.01; s.N=T/s.h; s.mc=10; t = [0:s.N]*s.h; 
s.m1=1; s.L1=1;    s.ell1=s.L1; s.I1=s.m1*s.ell1^2/3;
s.m2=0.5; s.L2=0.5;  s.ell2=s.L2; s.I2=s.m2*s.ell2^2/3; alpha=0.1;
s.B=[0; 0; 0; 1; 0; 0]; s.Q=eye(6); s.R = alpha^2; s.QT=diag([5 40 10 .1 60 10]);
s.x0=[0; pi; pi; 0; 0; 0]; g = 9.8;

%% Step 1: Optimizing U(t)
%Computing x_k with looped u_k
[u_k,x_k] = dual_inverted_pendulum(T,s);
state(:,1) = x_k(1:6,end);

% total_runs = 3;
% for i = 1:total_runs-1
%     [u_k(:,i+1),x_k] = dual_inverted_pendulum(T,u_k(:,end));
%     state(:,i+1) = x_k(1:6,end);
% end

% Determining A,N,E for each x_k
A = zeros(6,6,size(u_k,1));
X = zeros(6,6,size(u_k,1));
X0 = eye(6);
E = zeros(6,6,size(u_k,1));

for i = 1:length(x_k)
    A(:,:,i) = Compute_A(x_k(:,i),s);
    E(:,:,i) = Compute_E(x_k(:,i),s);
end

%% Step 2: Find K(t)
%Solve DRE using RK4 while marching backwards in time

X = RK4_backwards(X0,E,A,s.B,s.R,s.Q);
for i=1:size(X,3)
    K(:,:,i) = -s.R^(-1)*s.B'*X(:,:,i);
end

%% Step 3: Find L(t)
P = zeros(6,6,size(u_k,1));
P0 = eye(6);
C = diag([1 1 1 0 0 0]);
Q1 = eye(6);
Q2 = eye(6);

%Solve DRE using RK4 while marching forwards in time
P = RK4_forwards(P0,E,A,C,Q1,Q2);
for i = 1:size(P,3)
    L(:,:,i) = -P(:,:,i)*C'*Q2^(-1);
end

%% Step 4: Find K
E4 = [eye(3),zeros(3);...
      zeros(3),[s.mc+s.m1+s.m2,-s.m1*s.L1,-s.m2*s.L2;...
                -s.m1*s.L1,s.I1+s.m1*s.L1^2,0;...
                -s.m2*s.L2,0,s.I2+s.m2*s.L2^2]];
A4 = [zeros(3),eye(3);...
      [0,0,0;0,s.m1*g*s.L1,0;0,0,s.m2*g*s.L2],zeros(3)];

%Use ICARE
[X4,K4,L4] = icare(E4\A4,E4\s.B,s.Q,s.R);
K_new = -s.R^(-1)*s.B'*X4;
X_new = RK4_step4(x_k(1:6,end),E4,A4,s.B,K_new);

figure(1);
for i=1:6
    plot(linspace(0,1000,1001),X_new(i,:))
    hold on
end
title('LQR')
ylabel('X')
xlabel('t')
hold off

%% Step 5: Find L
[P,x,y] = icare(A4',C,Q1,Q2);
L_new = -P*C'*Q2^(-1);
x_hat = RK4_step5(x_k(1:6,end),A4,s.B,C,L_new,K_new,X_new);

figure(2);
for i=1:6
    plot(linspace(0,1000,1001),x_hat(i,:))
    hold on
end
title('LQG')
ylabel('X_hat')
xlabel('t')