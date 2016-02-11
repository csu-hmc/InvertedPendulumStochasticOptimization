function [K,Kd] = LQRcontroller(theta,params)

m = params.m;
l = params.l;
g = params.g;
J = m*l^2;

A = [0 1; m*g*l*sin(theta)/J 0];
B = [0;1/J];
R = 1;
Q = 1*eye(params.nstates);

Ks = lqr(A,B,Q,R);
K = Ks(1);
Kd = Ks(2);