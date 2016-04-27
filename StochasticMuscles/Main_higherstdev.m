rng('shuffle')

clear all
close all
clc

%Main program for pendulum swingup

% Declare parameters
params.m = 2;
params.l = 0.6;
params.g = 9.81;
params.T = 3; %Duration of motion
params.N = 300;
params.Nsamples = 1;
params.h = params.T/(params.N-1);

params.ndof = 1;
params.nmus = 2;
params.nstates = params.ndof*2+params.nmus*2;
params.ncontrols= 1;%params.nmus;
params.targetangle = pi/2;

params = getParams(params);
params.muscleparam = Muscleparams();
params.solver = 'IPOPT';
params.omega = 0*randn(params.ndof*2,params.N); %Added noise
params.W = 0;%1e-5;
params.method = 'BE';
params.epsilon = 1e-2;
params.xneutral = findNeutralstate(params);


% [X0, L, U] = getIniConBound(params);
% [~, params] = conjacstructure(L, U, params);

%load result 1 to 4
load('030516TempResults.mat');

params = result4.params;
randvals = params.omega/(1.2/sqrt(params.h));
stdev = 1.6/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result4);
result5 = Optimize(X0, L, U, params);

save('030716TempResults.mat', 'result1', 'result2', 'result3', 'result4', 'result5')

stdev = 2/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result5);
result6 = Optimize(X0, L, U, params);

save('030716TempResults.mat', 'result1', 'result2', 'result3', 'result4', 'result5', 'result6')

stdev = 2.5/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result6);
result7 = Optimize(X0, L, U, params);

save('030716TempResults.mat', 'result1', 'result2', 'result3', 'result4', 'result5', 'result6', 'result7')

stdev = 3/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result7);
result8 = Optimize(X0, L, U, params);

save('030716TempResults.mat', 'result1', 'result2', 'result3', 'result4', 'result5', 'result6', 'result7', 'result8')
