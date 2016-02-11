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
h = params.T/(params.NperSU+1);

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

% More noise
params.omega = 0.025*sqrt(h)*randn(params.nstates,params.N); 
[X0, L, U] = getIniConBound(params, result1);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result2';
params.warmstart = 1;
result2 = Optimize(X0, L, U, params);

% More noise
params.omega = 0.25*sqrt(h)*randn(params.nstates,params.N); 
[X0, L, U] = getIniConBound(params, result1);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result3';
params.warmstart = 1;
result3 = Optimize(X0, L, U, params);

% More noise
params.omega = 2.5*sqrt(h)*randn(params.nstates,params.N); 
[X0, L, U] = getIniConBound(params, result1);
[~, params] = conjacstructure(L, U, params);
params.snoptname = 'result4';
params.warmstart = 1;
result4 = Optimize(X0, L, U, params);

result(1) = result1;
result(2) = result2;
result(3) = result3;
result(4) = result4;
% result(5) = result5;
% result(6) = result6;
% result(7) = result7;