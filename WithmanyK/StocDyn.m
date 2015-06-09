function [f, dfdx, dfdxdot, dfdu, dfdK, dfdKd] = StocDyn(x, xdot, u, omega, params, K, Kd)

m = params.m;
l = params.l;
g = params.g;
torque = u;

theta = x(1);
dtheta = x(2);

J = m*l^2;
G = m*g*l*cos(theta);

a_x = [dtheta; -G/J];
Bu = [0; (torque+K*theta+Kd*dtheta)/J];
C = zeros(params.nstates);
C(2,2) = 1;

f = a_x + Bu + C*omega - xdot(1:params.nstates);
dfdx = [0 1;
    (m*g*l*sin(theta)+K)/J Kd/J];
dfdu = [0; 1/J];
dfdxdot = -1*eye(params.nstates);
dfdK = [0; theta/J];
dfdKd = [0; dtheta/J];