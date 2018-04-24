function [f, dfdx, dfdxdot, dfdu] = StocDyn(x, xdot, u, omega, params)

m = params.m;
l = params.l;
g = params.g;
nmus = params.nmus;
ndof = params.ndof;

theta = x(1);
dtheta = x(2);
b = 0;

%Get muscle information
[Fsee,dFseedx,fmus,dfmusdx,dfmusdxdot,dfmusdu] = getMusDyns(x,xdot,u,params);

% Find torque
d = params.muscleparam.d;
T = d*Fsee;
dTdx = d*dFseedx;

J = m*l^2;
G = m*g*l*cos(theta);

a_x = [dtheta; -G/J];
Bu = [0; (T-b*dtheta)/J];
C = zeros(params.ndof*2);
C(2,2) = 1;

f = a_x + Bu + C*omega - xdot(1:params.ndof*2);
f = [f;fmus];
dfdx = [[0 1; m*g*l*sin(theta)/J -b/J] zeros(params.ndof*2,params.nmus*2);dfmusdx];
dfdx(ndof*2,:) = dfdx(ndof*2,:)+dTdx/J;
dfdu = [zeros(ndof*2,params.nmus);dfmusdu];
dfdxdot = [-1*eye(params.ndof*2) zeros(params.ndof*2,params.nmus*2);dfmusdxdot];

% % scaling
% f(2:4) = f(2:4)/1e2;
% dfdx(2:4,:) = dfdx(2:4,:)/1e2;
% dfdu(2:4,:) = dfdu(2:4,:)/1e2;
% dfdxdot(2:4,:) = dfdxdot(2:4,:)/1e2;