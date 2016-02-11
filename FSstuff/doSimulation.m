function [x1,u1,K1,timesFS,params] = doSimulation(params, X, stdev, feedback)

params.Ntimes = params.NperSU;
params.feedback = feedback;
xu = reshape(X(1:end-params.Ks),params.nvarpernode, params.N);

xperSU = zeros(params.nstates, params.NperSU, params.NSU);
uperSU = zeros(params.ncontrols,params.NperSU, params.NSU);
for i = 1:params.NSU
    xperSU(:,:,i) = xu(1:params.nstates,(i-1)*params.NperSU+(1:params.NperSU));
    uperSU(:,:,i) = xu(params.nstates+1:end,(i-1)*params.NperSU+(1:params.NperSU));
end
meanx = mean(xperSU,3);
meanu = mean(uperSU,3);
fu = meanu;
K = X(end-1);
Kd= X(end);
u0only = fu+K*meanx(1,:)+Kd*meanx(2,:);

if params.Ks == 120
    K = X(end-119:end-60);
    Kd = X(end-59:end);
    u0only = fu+K'.*meanx(1,:)+Kd'.*meanx(2,:);
end

timesDC = linspace(0,10,params.NperSU);
timesFS = linspace(0,10,params.Ntimes);
for i = 1:params.nstates
    params.xdes(i,:) = interp1(timesDC, meanx(i,:), timesFS, 'spline');
end
params.u0 = interp1(timesDC, u0only, timesFS, 'spline');
params.omega = stdev*randn(params.nstates,params.Ntimes);
[x1,u1,K1,S1] = ForDynLQR(params);