rng('shuffle')

clear all
close all
clc

%Main program for pendulum swingup

params.asat = 0;

% Declare parameters
params.m = 2;
params.l = 0.6;
params.g = 9.81;
params.T = 5; %Duration of motion
params.N = 100;
params.Nsamples = 1;
params.h = params.T/(params.N-1);

params.ndof = 1;
params.nmus = 2;
params.nstates = params.ndof*2+params.nmus*2;
params.ncontrols= 2;%params.nmus;
params.Ks = params.ncontrols*2;
params.targetangle = pi/2;

params = getParams(params);
params.muscleparam = Muscleparams();
params.solver = 'IPOPT';
params.omega = 0*randn(params.ndof*2,params.N); %Added noise
params.W = 0;%1e-5;
params.method = 'BE';
params.epsilon = 1e-1;
params.xneutral = findNeutralstate(params);

[X0, L, U] = getIniConBound(params);
[~, params] = getConjacstructure(L, U, params);
% 
% [grad,grad_num,cjac, cjac_num] = checkDerivatives(X0, params);%+randn(size(X0))
% keyboard;

%First result witout noise for desired trajectory
result1 = Optimize(X0, L, U, params);

% keyboard
% plotresult(result1);
% pause(1e-6)
% [X0, L, U] = getIniConBound(params,1,result1);
[~, params] = getConjacstructure(L, U, params);
% 
% %First result witout noise for desired trajectory
% result1n = Optimize(X0, L, U, params);

params.Nsamples = 15;
params = getParams(params);

Ws = [1e-2 1e-1 1 10 100];
% for i = 1:length(Ws)
%     params.W = Ws(i);
    % Now add noise and increase time
    stdev = 0.01/sqrt(params.h);
    randvals = randn(params.ndof*2,params.N*params.Nsamples); %Added noise
    params.omega = stdev*randvals;

    [X0, L, U] = getIniConBound(params,0,result1);%result1);
    [~, params] = getConjacstructure(L, U, params);
    
    result2 = Optimize(X0, L, U, params);

    keyboard
    [X0, L, U] = getIniConBound(params,1,result2);
    result2n = Optimize(X0, L, U, params);

    stdev = 0.05/sqrt(params.h);
    params.omega = stdev*randvals;

    [X0, L, U] = getIniConBound(params,0,result2);
    result3 = Optimize(X0, L, U, params);

    keyboard
%     [X0, L, U] = getIniConBound(params,1,result2n);
%     result3n = Optimize(X0, L, U, params);

    stdev = 0.1/sqrt(params.h);
    params.omega = stdev*randvals;

    [X0, L, U] = getIniConBound(params,0,result3);
    result4 = Optimize(X0, L, U, params);

%     keyboard
%     [X0, L, U] = getIniConBound(params,1,result3n);
%     result4n = Optimize(X0, L, U, params);

    stdev = 0.5/sqrt(params.h);
    params.omega = stdev*randvals;
    
    [X0, L, U] = getIniConBound(params,0,result4);
    result5 = Optimize(X0, L, U, params);
   
%     [X0, L, U] = getIniConBound(params,1,result4n);
%     result5n = Optimize(X0, L, U, params);

    stdev = 1/sqrt(params.h);
    params.omega = stdev*randvals;
    
    [X0, L, U] = getIniConBound(params,0,result5);
    result6 = Optimize(X0, L, U, params);
   
%     [X0, L, U] = getIniConBound(params,1,result5n);
%     result6n = Optimize(X0, L, U, params);
    
%     filename = ['042516NoSatur' num2str(i) '.mat'];
%     save(filename, 'result1', 'result2', 'result3', 'result4', 'result5', 'result2n', 'result3n', 'result4n', 'result5n', 'result6', 'result6n')
% % 
% end
