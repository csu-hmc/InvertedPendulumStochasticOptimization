function [x,u] = testTrajs(stdev, params,result)

%Now find feedback gains
X = result.X;
xu = reshape(X(1:end-2),params.nvarpernode, params.N);

xperSU = zeros(params.nstates, params.NperSU, params.NSU);
uperSU = zeros(params.ncontrols,params.NperSU, params.NSU);
for i = 1:params.NSU
    xperSU(:,:,i) = xu(1:params.nstates,(i-1)*params.NperSU+(1:params.NperSU));
    uperSU(:,:,i) = xu(params.nstates+1:end,(i-1)*params.NperSU+(1:params.NperSU));
end
meanx = mean(xperSU,3);
meanu = mean(uperSU,3);

timesDC = linspace(0,10,params.NperSU);
timesFS = linspace(0,10,params.Ntimes);
for i = 1:params.nstates
    params.xdes(i,:) = interp1(timesDC, meanx(i,:), timesFS, 'spline');
end
params.u0 = interp1(timesDC, meanu, timesFS, 'spline');
params.omega = stdev*randn(params.nstates,params.Ntimes);

[x,u] = ForDynLQR(params);