function [f, dfdx, dfdxdot, dfdu] = StocDyn(x, xdot, u, omega, params)

m = params.m;
l = params.l;
g = params.g;
torque = u;
Omega_y = params.Omega_y;
Omega_m = params.Omega_m;

theta = x(1);
dtheta = x(2);
mhat = x(3);
Sigma = x(4);

J = m*l^2;
G = m*g*l*cos(theta);

a_x = [dtheta; -G/J; 0; Omega_m-Sigma^2*theta^2/Omega_y];
Bu = [0; torque/J; 0; 0];
C = zeros(4);
C(3,3) = Sigma*theta/Omega_y;

f = a_x + Bu + C*omega - xdot;
dfdx = [0 1 0 0;
    m*g*l*sin(theta)/J 0 0 0;
    Sigma*omega(3)/Omega_y 0 0 theta*omega(3)/Omega_y;
    -2*Sigma^2*theta/Omega_y 0 0 -2*Sigma*theta^2/Omega_y];
dfdu = [0; 1/J; 0 ; 0];
dfdxdot = -1*eye(4);