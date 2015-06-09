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

params = getParams(params);
params.omega = zeros(params.nstates,params.N);

[X0, L, U] = getIniConBound(params);

[~, params] = conjacstructure(L, U, params);

%First result witout noise for desired trajectory
result1 = Optimize(X0, L, U, params);
params1 = params;

% Now add more swingups and noise
params.NSU = 100;
params = getParams(params);

params.omega = 0.001*randn(params.nstates,params.N); %Added noise, low pass filtered
% [b,a] = butter(2,0.1);
% params.omega = filter(b,a,omega);

[X0, L, U] = getIniConBound(params, result1);
[~, params] = conjacstructure(L, U, params);

result2 = Optimize(X0, L, U, params);
params2 = params;
% plotresult(result, params);

% More noise
X0 = result2.X;
params.omega = 0.01*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result3 = Optimize(X0, L, U, params);
params3 = params;

% profile on
% More noise
X0 = result3.X;
params.omega = 0.05*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result4 = Optimize(X0, L, U, params);
params4 = params;
% profile viewer

% More noise
X0 = result4.X;
params.omega = 0.1*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result5 = Optimize(X0, L, U, params);
params5 = params;

% More noise
X0 = result5.X;
params.omega = 0.3*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result6 = Optimize(X0, L, U, params);
params6 = params;

% More noise
X0 = result6.X;
params.omega = 1*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result7 = Optimize(X0, L, U, params);
params7 = params;