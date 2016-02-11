function [x,uout,K,S] = ForDynLQR_td(params)

N = params.Ntimes;
T = params.T;
h = T/(N-1);
u0 = params.u0;
xdes = params.xdes;

x = zeros(params.nstates,N);
uout = zeros(params.ncontrols,N);
u = u0;
x(:,1) = [-pi/2;0];
uout(1) = u0(1);

omega = params.omega;

R = params.R;%1;
Q = params.Q;%10*eye(params.nstates);
[K,S,E] = dplqr(Q, R, xdes, N, params);
if params.feedback ~= 1
    K = zeros(size(K));
end
uout(1) = u0(1);

% Newton Step for midpoint Euler
maxIters = 100;
tolerance = 1e-4;

% No feedback yet in first time step
x(:,2) = x(:,1);
tol = 1;
iters = 0;
while and(tol>tolerance,iters<maxIters)
    x1 = x(:,1);
    x2 = x(:,2);
    u1 = u(1);
    u2 = u(2)-K(:,:,2)*[x1(1)-xdes(1,1);x1(2)-xdes(2,1)];
    [f, dfdx, dfdxdot, dfdu] = StocDynFS((x1+x2)/2, (x2-x1)/h, (u1+u2)/2, omega(:,1), params);
    dfdx2 = 1/2*dfdx+1/h*dfdxdot+1/2*dfdu*-K(:,:,2);
    xold = x(:,2);
    x(:,2) = x(:,2)-dfdx2\f;
    iters = iters+1;
    tol = max(abs(xold-x(:,2)));
end

for i = 3:N
    x(:,i) = x(:,i-1);
    tol = 1;
    iters = 0;
    while and(tol>tolerance,iters<maxIters)
        x1 = x(:,i-1);
        x2 = x(:,i);
        u1 = u(i-1)-K(:,:,i-1)*[x1(1)-xdes(1,i-2);x1(2)-xdes(2,i-2)];
        u2 = u(i)-K(:,:,i)*[x2(1)-xdes(1,i-1);x2(2)-xdes(2,i-1)];
        [f, dfdx, dfdxdot, dfdu] = StocDynFS((x1+x2)/2, (x2-x1)/h, (u1+u2)/2, omega(:,i-1), params);
        dfdx2 = 1/2*dfdx+1/h*dfdxdot+1/2*dfdu*-K(:,:,i);
        xold = x(:,i);
        x(:,i) = x(:,i)-dfdx2\f;
        iters = iters+1;
        tol = max(abs(xold-x(:,i)));
    end
    uout(i) = u2;
end
    

params.omega = omega;