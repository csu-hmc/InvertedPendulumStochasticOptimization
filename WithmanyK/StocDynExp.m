function xdot = StocDynExp(t, x, u, params)

m = params.m;
l = params.l;
g = params.g;

theta = x(1);
dtheta = x(2);

J = m*l^2;
G = m*g*l*cos(theta);

a_x = [dtheta; -G/J];
Bu = [0; u/J];
C = zeros(2);
C(2,2) = 1;

omega = 0*randn(2,1);
xdot = a_x + Bu + C*omega;