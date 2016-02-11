function [x,uout] = ForDynKs(params)

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

%4th order Runge-Kutta
for i = 2:N
    params.omega = omega(:,i);
    [k1,uout(i)] = DynK(x(:,i-1),u(i-1),xdes(:,i-1),params);
    k2 = DynK(x(:,i-1)+h/2*k1,u(i-1),xdes(:,i-1),params);
    k3 = DynK(x(:,i-1)+h/2*k2,u(i-1),xdes(:,i-1),params);
    k4 = DynK(x(:,i-1)+h*k3,u(i-1),xdes(:,i-1),params);
    x(:,i) = x(:,i-1)+h/6*(k1+2*k2+2*k3+k4);
end

params.omega = omega;