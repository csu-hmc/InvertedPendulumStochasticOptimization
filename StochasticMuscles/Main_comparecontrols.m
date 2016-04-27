rng('shuffle')

% clear all
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
params.W = 000;%1e-5;
params.method = 'BE';
params.epsilon = 1e-2;
params.xneutral = findNeutralstate(params);


[X0, L, U] = getIniConBound(params);
[~, params] = conjacstructure(L, U, params);
% 
% [grad,grad_num,cjac, cjac_num] = checkDerivatives(X0+randn(size(X0)), params);
% keyboard;

%First result witout noise for desired trajectory
result1 = Optimize(X0, L, U, params);

keyboard
% plotresult(result1);
% pause(1e-6)
% [X0, L, U] = getIniConBound(params,1,result1);
% [~, params] = conjacstructure(L, U, params);
% 
% %First result witout noise for desired trajectory
% result1n = Optimize(X0, L, U, params);

params.Nsamples = 5;
params = getParams(params);
params.Jnnz = params.Jnnz*params.Nsamples;
% for i = 1:10
    % Now add noise and increase time
    stdev = .5/sqrt(params.h);
    randvals = randn(params.ndof*2,params.N*params.Nsamples); %Added noise
    params.omega = stdev*randvals;

    [X0, L, U] = getIniConBound(params,0,result1);%result1);
%     [~, params] = conjacstructure(L, U, params);
    
%     [grad,grad_num,cjac, cjac_num] = checkDerivatives(X0+randn(size(X0)), params);
%     keyboard;

    result2 = Optimize(X0, L, U, params);
    
%     save('032516highmass.mat', 'result1', 'result2')

    keyboard;
%     [X0, L, U] = getIniConBound(params,1,result2);
%     result2n = Optimize(X0, L, U, params);

    stdev = 1/sqrt(params.h);
    params.omega = stdev*randvals;

    [X0, L, U] = getIniConBound(params,0,result2);
    result3 = Optimize(X0, L, U, params);

%     save('032516highmass.mat', 'result1', 'result2', 'result3')
    keyboard;
%     [X0, L, U] = getIniConBound(params,1,result3);
%     result3n = Optimize(X0, L, U, params);

    stdev = 1.5/sqrt(params.h);
    params.omega = stdev*randvals;

    [X0, L, U] = getIniConBound(params,0,result3);
    result4 = Optimize(X0, L, U, params);

%     save('032516highmass.mat', 'result1', 'result2', 'result3', 'result4')
%     [X0, L, U] = getIniConBound(params,1,result4);
%     result4n = Optimize(X0, L, U, params);
keyboard
    stdev = 2/sqrt(params.h);
    params.omega = stdev*randvals;
    
    [X0, L, U] = getIniConBound(params,0,result4);
    result5 = Optimize(X0, L, U, params);

% save('032516highmass.mat', 'result1', 'result2', 'result3', 'result4', 'result5')


    stdev = 3/sqrt(params.h);
    params.omega = stdev*randvals;
    
    [X0, L, U] = getIniConBound(params,0,result5);
    result6 = Optimize(X0, L, U, params);

% save('032516highmass.mat', 'result1', 'result2', 'result3', 'result4', 'result5', 'result6')

% end

    stdev = 10/sqrt(params.h);
    params.omega = stdev*randvals;
    
    [X0, L, U] = getIniConBound(params,0,result6);
    result7 = Optimize(X0, L, U, params);

% save('032516highmass.mat', 'result1', 'result2', 'result3', 'result4', 'result5', 'result6', 'result7')

% end
