function [xdot,uout] = DynK(x,u,xdes,params)

m = params.m;
l = params.l;
g = params.g;
torque = u;

theta = x(1);
dtheta = x(2);

J = m*l^2;
G = m*g*l*cos(theta);
K = params.K;
Kd = params.Kd;

uout = torque-K*(theta-xdes(1))-Kd*(dtheta-xdes(2));
a_x = [dtheta; -G/J];
Bu = [0; uout/J];
C = zeros(params.nstates);
C(2,2) = 1;

xdot = a_x + Bu + C*params.omega;