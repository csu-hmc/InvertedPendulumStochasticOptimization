clear all
close all
clc

rng('shuffle')

% Main program for pendulum swingup
% global params
% Declare parameters
params.m = 5;
params.l = 1.2;
params.g = 9.81;
params.T = 10; %Duration of motion

params.nstates  = 2;
params.ncontrols= 1;

params.NperSU = 60;
params.NSU = 1; %number of swingups
params.Ks = params.NperSU*2;
params.ineq = 0;
params.optimizer = 'IPOPT';
params.acc = 0;
params.obj = 0;

params = getParams(params);
params.omega = zeros(params.nstates,params.N);

[X0, L, U] = getIniConBound(params);

[J, params] = conjacstructure(L, U, params);
% [grad,grad_num,cjac, cjac_num] = checkDerivatives(params,L,U);
% keyboard

%First result witout noise for desired trajectory
params.snoptname = 'result1';
params.warmstart = 1;
result1 = Optimize(X0, L, U, params);
disp('First one done')
result1

% Now add more swingups and noise
params.NSU = 50;
params = getParams(params);

params.omega = 0.001*randn(params.nstates,params.N); %Added noise, low pass filtered 0.001*

[X0, L, U] = getIniConBound(params, result1);
[~, params] = conjacstructure(L, U, params);

% [grad,grad_num,cjac, cjac_num] = checkDerivatives(params,X0); %+10*randn(size(L)));
% keyboard

params.snoptname = 'result2';
params.warmstart = 1;
result2 = Optimize(X0, L, U, params);
disp('Second one done')
result2


% % debuggingcode

% More noise
X0 = result2.X;
params.omega = 0.01*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result2);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result3';
params.warmstart = 1;
result3 = Optimize(X0, L, U, params);
disp('Third one done')
result3
% 
% More noise
X0 = result3.X;
params.omega = 0.1*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result3);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result4';
params.warmstart = 1;
result4 = Optimize(X0, L, U, params);
disp('Fourth one done')
result4
% % 
% More noise
X0 = result4.X;
params.omega = 0.5*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result4);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result5';
params.warmstart = 1;
result5 = Optimize(X0, L, U, params);
disp('Fifth one done')
result5

% More noise
X0 = result5.X;
params.omega = 1*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result5);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result6';
params.warmstart = 1;
result6 = Optimize(X0, L, U, params);
disp('Sixth one done')
result6

% More noise
X0 = result6.X;
params.omega = 1.2*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result6);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result7';
params.warmstart = 1;
result7 = Optimize(X0, L, U, params);
disp('Seventh one done')
result7
plotresult(result7);

result(1) = result1;
result(2) = result2;
result(3) = result3;
result(4) = result4;
result(5) = result5;
result(6) = result6;
result(7) = result7;
for i = 1:length(result)
    obj(i) = result(i).obj;
    [xnow, u0,K,Kd] = plotresult(result(i), 0);
        for k = 1:length(u0)
        uall(k) = findTorque(u0(k), [K(k);Kd(k)], xnow(:,k));
        end
    meanx(i,:) = xnow(1,:);
    rmsu = rms(uall);
    RMSu(i) = mean(rmsu(:));
end
figure
plot(meanx')