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
params.T = 20; %Duration of motion
% params.K = -2; %Feedback gain 

params.nstates  = 2;
params.optstates= 4;
params.ncontrols= 1;

params.NperSU = 60;
params.NSU = 1; %number of swingups

params = getParams(params);
params.omega = zeros(params.nstates,params.N);

[X0, L, U] = getIniConBound(params);

[~, params] = conjacstructure(L, U, params);
% [grad,grad_num,cjac, cjac_num] = checkDerivatives(L, U, params);

%First result witout noise for desired trajectory
result1 = Optimize(X0, L, U, params);
params1 = params;

% Now add more swingups and noise
params.NSU = 15;
params = getParams(params);

params.omega = 0.5*randn(params.nstates,params.N); %Added noise, low pass filtered
% [b,a] = butter(2,0.1);
% params.omega = filter(b,a,omega);

[X0, L, U] = getIniConBound(params, result1);
[~, params] = conjacstructure(L, U, params);

result = Optimize(X0, L, U, params);
plotresult(result, params);