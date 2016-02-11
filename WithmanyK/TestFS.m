clear all
close all
clc

load('ResultUgh.mat')
params = result2.params;
stdev = 0;%0.001;
X = result2.X;

X1 = X(1:params.nvarSU1);
x1 = reshape(X1,params.nvarpernode1,params.NperSU);
xperSU(:,:,1) = x1(1:2,:);
u0 = x1(3,:);
K = x1(4,:);
Kd = x1(5,:);
if params.NSU > 1
    Xelse = X(params.nvarSU1+1:end);
    xelse = reshape(Xelse,params.nvarpernode,params.NperSU,params.NSU-1);
    xperSU(:,:,2:params.NSU) = xelse;
end

params.Ntimes = params.NperSU;
omega = stdev*randn(params.nstates,params.Ntimes);
N = params.Ntimes;
T = params.T;
h = T/(N-1);

x = zeros(params.nstates,N);
uout = zeros(params.ncontrols,N);
x(:,1) = [-pi/2;0];

% Using time dependent gains from optimization
uout(1) = u0(1)+[K(1) Kd(1)]*[x(1,1);x(2,1)];

% Newton Step for midpoint Euler
maxIters = 100;
tolerance = 1e-4;
for i = 2:N
    x(:,i) = x(:,i-1);
    tol = 1;
    iters = 0;
    while and(tol>tolerance,iters<maxIters)
        x1 = x(:,i-1);
        x2 = x(:,i);
                 
        u1 = u0(i-1)+[K(i-1) Kd(i-1)]*[x1(1);x1(2)];
        u2 = u0(i)+[K(i) Kd(i)]*[x2(1);x2(2)];
        [f, dfdx, dfdxdot, dfdu] = StocDynFS((x1+x2)/2, (x2-x1)/h, (u1+u2)/2, omega(:,i-1), params);
        dfdx2 = 1/2*dfdx+1/h*dfdxdot+1/2*dfdu*[K(i) Kd(i)];

        xold = x(:,i);
        x(:,i) = x(:,i)-dfdx2\f;
        iters = iters+1;
        tol = max(abs(xold-x(:,i)));
    end
    uout(i) = u2;
end   

