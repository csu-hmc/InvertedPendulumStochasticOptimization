clear all
close all
clc

%Main program for pendulum swingup

% Declare parameters
params.m = 5;
params.l = 1.2;
params.g = 9.81;
params.ktheta = 0.01;
params.kthetadot = 0.01;
params.T = 10; %Duration of motion

params.nstates  = 2;
params.ncontrols= 1;

params.NperSU = 60;
params.NSU = 1; %number of swingups
params.ineq = 0; %Enable inequality constraints

params = getParams_Ktest(params);
params.omega = zeros(params.nstates,params.N);
params.solver = 'IPOPT';

[X0, L, U] = getIniConBound(params);

[~, params] = conjacstructure(L, U, params);
% [grad,grad_num,cjac, cjac_num] = checkDerivatives(X0, params);

%First result witout noise for desired trajectory
result1 = Optimize(X0, L, U, params);

% Now add more swingups and noise
params.NSU = 30;
params.ineq = 0;
params = getParams_Ktest(params);

params.omega = 0.001*randn(params.nstates,params.N); %Added noise

[X0, L, U] = getIniConBound(params, result1);
[~, params] = conjacstructure(L, U, params);
% [grad,grad_num,cjac, cjac_num] = checkDerivatives(X0, params);
% keyboard;

result2 = Optimize(X0, L, U, params);
% plotresult(result2, 1);
% 
% More noise
X0 = result2.X;
params.omega = 0.01*randn(params.nstates,params.N); %Added noise
[~, params] = conjacstructure(L, U, params);
result3 = Optimize(X0, L, U, params);
% 
% More noise
X0 = result3.X;
params.omega = 0.1*randn(params.nstates,params.N); %Added noise
[~, params] = conjacstructure(L, U, params);
result4 = Optimize(X0, L, U, params);

% More noise
X0 = result4.X;
params.omega = 0.5*randn(params.nstates,params.N); %Added noise
[~, params] = conjacstructure(L, U, params);
result5 = Optimize(X0, L, U, params);

% More noise
X0 = result5.X;
params.omega = 0.7*randn(params.nstates,params.N); %Added noise
[~, params] = conjacstructure(L, U, params);
result6 = Optimize(X0, L, U, params);

% More noise
X0 = result6.X;
params.omega = 1*randn(params.nstates,params.N); %Added noise
[~, params] = conjacstructure(L, U, params);
result7 = Optimize(X0, L, U, params);