function [x,uout,K,S] = ForDynLQR(params)

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
% % Using LQR
% uout(1) = u0(1)-K(:,:,1)*[x(1,1)-xdes(1,1);x(2,1)-xdes(2,1)];

% Using time dependent gains from optimization
uout(1) = u0(1)+[params.K(1) params.Kd(1)]*[x(1,1)-xdes(1,1);x(2,1)-xdes(2,1)];

% % Using the two gains from optimization
% uout(1) = u0(1)+[params.K params.Kd]*[x(1,1)-xdes(1,1);x(2,1)-xdes(2,1)];
% K = [params.K params.Kd];
% S= 0;

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
        
        % Using LQR
        u1 = u(i-1)-K(:,:,i-1)*[x1(1)-xdes(1,i-1);x1(2)-xdes(2,i-1)];
        u2 = u(i)-K(:,:,i)*[x2(1)-xdes(1,i);x2(2)-xdes(2,i)];
                 
%         % Using optimized gains
%         if params.Ks == 120
%             u1 = u(i-1)+[params.K(i-1) params.Kd(i-1)]*[x1(1)-xdes(1,i-1);x1(2)-xdes(2,i-1)];
%             u2 = u(i)+[params.K(i) params.Kd(i)]*[x2(1)-xdes(1,i);x2(2)-xdes(2,i)];
%         else            
%             u1 = u(i-1)+[params.K params.Kd]*[x1(1)-xdes(1,i-1);x1(2)-xdes(2,i-1)];
%             u2 = u(i)+[params.K params.Kd]*[x2(1)-xdes(1,i);x2(2)-xdes(2,i)];
%         end
        [f, dfdx, dfdxdot, dfdu] = StocDynFS((x1+x2)/2, (x2-x1)/h, (u1+u2)/2, omega(:,i-1), params);
        dfdx2 = 1/2*dfdx+1/h*dfdxdot+1/2*dfdu*-K(:,:,i);
%         if params.Ks == 120
%             dfdx2 = 1/2*dfdx+1/h*dfdxdot+1/2*dfdu*[params.K(i) params.Kd(i)];
%         else
%             dfdx2 = 1/2*dfdx+1/h*dfdxdot+1/2*dfdu*[params.K params.Kd];
%         end

        xold = x(:,i);
        x(:,i) = x(:,i)-dfdx2\f;
        iters = iters+1;
        tol = max(abs(xold-x(:,i)));
    end
    uout(i) = u2;
end   

params.omega = omega;


% %4th order Runge-Kutta
% for i = 2:N
%     params.omega = omega(:,i-1);
%     params.K = K(:,:,i);
%     [k1,uout(i)] = DynLQR(x(:,i-1),u(i-1),xdes(:,i-1),params);
%     k2 = DynLQR(x(:,i-1)+h/2*k1,u(i-1),xdes(:,i-1),params);
%     k3 = DynLQR(x(:,i-1)+h/2*k2,u(i-1),xdes(:,i-1),params);
%     k4 = DynLQR(x(:,i-1)+h*k3,u(i-1),xdes(:,i-1),params);
%     x(:,i) = x(:,i-1)+h/6*(k1+2*k2+2*k3+k4);
% end
