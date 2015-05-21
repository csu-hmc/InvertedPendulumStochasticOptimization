function xdot = StocDyn(t, x, u, omega, params)

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

xdot = a_x + Bu + C*omega;