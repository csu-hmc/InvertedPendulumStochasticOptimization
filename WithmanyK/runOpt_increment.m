function result = runOpt(NSU, optimizer)

% if nargin == 1
%     optimizer = 'SNOPT';
% end
% 
% global params
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
params.optimizer = optimizer;
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

% Now add more swingups and noise
params.NSU = NSU;
params = getParams(params);

params.omega = 0.001*randn(params.nstates,params.N); %Added noise, low pass filtered

[X0, L, U] = getIniConBound(params, result1);
[~, params] = conjacstructure(L, U, params);

params.snoptname = 'result2';
params.warmstart = 1;
result2 = Optimize(X0, L, U, params);

% More noise
X0 = result2.X;
params.omega = 0.01*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result2);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result3';
params.warmstart = 1;
result3 = Optimize(X0, L, U, params);

% More noise
X0 = result3.X;
params.omega = 0.1*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result3);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result4';
params.warmstart = 1;
result4 = Optimize(X0, L, U, params);

% More noise
X0 = result4.X;
params.omega = 0.3*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result4);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result5';
params.warmstart = 1;
result5 = Optimize(X0, L, U, params);

% More noise
X0 = result5.X;
params.omega = 0.5*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result5);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result6';
params.warmstart = 1;
result6 = Optimize(X0, L, U, params);

% More noise
X0 = result6.X;
params.omega = 1*randn(params.nstates,params.N); %Added noise, low pass filtered
[X0, L, U] = getIniConBound(params, result6);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result7';
params.warmstart = 1;
result7 = Optimize(X0, L, U, params);

result(1) = result1;
result(2) = result2;
result(3) = result3;
result(4) = result4;
result(5) = result5;
result(6) = result6;
result(7) = result7;