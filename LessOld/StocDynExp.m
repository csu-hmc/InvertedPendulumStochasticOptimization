function xdot = StocDynExp(t, x, params)

m = params.m;
l = params.l;
g = params.g;
torque = 1;% u(1);
torquedes = 1;%u(2);
K = params.K;

theta = x(1);
dtheta = x(2);
thetades = x(3);
dthetades = x(4);

J = m*l^2;
G = m*g*l*cos(theta);
Gdes = m*g*l*cos(thetades);

a_x = [dtheta; -G/J; dthetades; -Gdes/J];
Bu = [0; torque/J+K*(theta-thetades); 0; torquedes/J];
C = zeros(4);
C(2,2) = 1;

omega = 100*randn(4,1);
xdot = a_x + Bu + C*omega;