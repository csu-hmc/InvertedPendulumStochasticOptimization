function [f, Bu, dfdx, dfdxdot, dfdu, dfdK] = StocDyn(x, xdot, u, omega, params)

m = params.m;
l = params.l;
g = params.g;

theta = x(1);
dtheta = x(2);

J = m*l^2;
G = m*g*l*cos(theta);

a_x = [dtheta; -G/J];
dadx = [0 1; m*g*l*sin(theta)/J 0];
Bu = [0; u/J];
dBdu = [0; 1/J];
% dBdx = [0 0; 1/J*dudx];
% dBdK = [0 0; 1/J*dudK];
C = zeros(params.nstates);
C(2,2) = 1;

f = a_x + Bu + C*omega - xdot(1:params.nstates);
dfdx = dadx;%+dBdx;
dfdu = dBdu;
dfdxdot = -1*eye(params.nstates);
dfdK = zeros(2);%dBdK;