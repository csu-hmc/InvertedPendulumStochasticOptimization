clear all
close all
clc

%Main program to test swingup result

load('swingupresult.mat');

X = result4.X;
xx = reshape(X(1:end),params4.nvarpernode, params4.N);
xperSU = xx(1:params4.optstates,1:params4.NperSU);
uperSU = xx(params4.optstates+1:end,1:params4.NperSU);
timestep = params4.T/(params4.NperSU-1);

x = zeros(2,params4.NperSU);
t = zeros(params4.NperSU,1);

uu = uperSU(1)+sum(xperSU(1:2,1).*xperSU(3:4,1));
[T,Y] = ode45(@(t,x) StocDynExp(t, x, uu, params4),[0 timestep],xperSU(1:2,1));
t(1) = T(end);
x(:,1) = Y(end,:);
for i = 2:params4.NperSU
    xini = x(:,i-1);
    uu = uperSU(i)+sum(xini.*xperSU(3:4,i));
    [T,Y] = ode45(@(t,x) StocDynExp(t, x, uu, params4),[t(i-1) t(i-1)+timestep],xini);
    t(i) = T(end);
    x(:,i) = Y(end,:);
end

plot(t,x(1,:))
