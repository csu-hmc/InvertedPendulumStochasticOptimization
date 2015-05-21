function [f, dfdx, dfdxdot, dfdu, dfdK, dfdKd] = StocDyn(x, xdot, u, K, Kd, omega, params)

m = params.m;
l = params.l;
g = params.g;
torque = u(1);
torquedes = u(2);

theta = x(1);
dtheta = x(2);
thetades = x(3);
dthetades = x(4);

J = m*l^2;
G = m*g*l*cos(theta);
Gdes = m*g*l*cos(thetades);

a_x = [dtheta; -G/J; dthetades; -Gdes/J];
Bu = [0; (torque+K*(theta-thetades)+Kd*(dtheta-dthetades))/J; 0; torquedes/J];
C = zeros(4);
C(2,2) = 1;

f = a_x + Bu + C*omega - xdot;
dfdx = [0 1 0 0;
    (m*g*l*sin(theta)+K)/J Kd/J -K/J -Kd/J;
    0 0 0 1;
    0 0 m*g*l*sin(thetades)/J 0];
dfdu = [0 0; 1/J 0; 0 0; 0 1/J];
dfdxdot = -1*eye(4);
dfdK = [0; (theta-thetades)/J; 0; 0];
dfdKd = [0; (dtheta-dthetades)/J; 0; 0];