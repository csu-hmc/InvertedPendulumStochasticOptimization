rng('shuffle')

clear all
close all
clc

%Main program for pendulum swingup
% global params

% Declare parameters
params.m = 2;
params.l = 0.6;
params.g = 9.81;
params.T = 5; %Duration of motion
params.N = 300;
params.h = params.T/(params.N-1);

params.ndof = 1;
params.nmus = 2;
params.nstates = params.ndof*2+params.nmus*2;
params.ncontrols= params.nmus;
params.Ks = params.ncontrols*2;
params.targetangle = pi/2;

params.muscleparam = Muscleparams();
params.solver = 'IPOPT';
params.W = 0;%1e-5;
params.method = 'BE';
params.epsilon = 1e-2;
params.xneutral = findNeutralstate(params);
 
lengthi = i;

Ws = 100;%[0.01 0.05 0.1 0.5 1 5 10 50 100];
stdevs = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];% 1.5 2 5];

params.W = Ws;

%First result witout noise for desired trajectory
params.Nsamples = 1;
params = getParams(params);
params.omega = 0*randn(params.ndof*2,params.N*params.Nsamples); %Added noise
[X0, L, U] = getIniConBound(params,0);
[~, params] = getConjacstructure(L, U, params);
%     [grad,grad_num,cjac, cjac_num] = checkDerivatives(X0+randn(size(X0)), params);%
%     keyboard;
params.snoptname = 'result1';
params.warmstart = 0;
result1 = Optimize(X0, L, U, params);

params.warmstart = 1;

% Now add noise and increase time
stdev = 0.1/sqrt(params.h);
randvals = randn(params.ndof*2,params.N*params.Nsamples); %Added noise
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result1);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result1.zl;
    params.zu = result1.zu;
    params.lambda = result1.lambda;
end

% Now add noise and increase time
stdev = 0.2/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result2);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result2.zl;
    params.zu = result2.zu;
    params.lambda = result2.lambda;
end


result3 = Optimize(X0, L, U, params);

% Now add noise and increase time
stdev = 0.3/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result3);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result3.zl;
    params.zu = result3.zu;
    params.lambda = result3.lambda;
end

result4 = Optimize(X0, L, U, params);

% Now add noise and increase time
stdev = 0.4/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result4);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result4.zl;
    params.zu = result4.zu;
    params.lambda = result4.lambda;
end

result5 = Optimize(X0, L, U, params);

% Now add noise and increase time
stdev = 0.5/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result5);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result5.zl;
    params.zu = result5.zu;
    params.lambda = result5.lambda;
end

result6 = Optimize(X0, L, U, params);

% Now add noise and increase time
stdev = 0.6/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result6);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result6.zl;
    params.zu = result6.zu;
    params.lambda = result6.lambda;
end

result7 = Optimize(X0, L, U, params);

% Now add noise and increase time
stdev = 0.7/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result7);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result7.zl;
    params.zu = result7.zu;
    params.lambda = result7.lambda;
end

result8 = Optimize(X0, L, U, params);

% Now add noise and increase time
stdev = 0.8/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result8);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result8.zl;
    params.zu = result8.zu;
    params.lambda = result8.lambda;
end

result9 = Optimize(X0, L, U, params);

% Now add noise and increase time
stdev = 0.9/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result9);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result9.zl;
    params.zu = result9.zu;
    params.lambda = result9.lambda;
end

result10 = Optimize(X0, L, U, params);

% Now add noise and increase time
stdev = 1/sqrt(params.h);
params.omega = stdev*randvals;

[X0, L, U] = getIniConBound(params,0,result10);%result1);
[~, params] = getConjacstructure(L, U, params);

if strcmp(params.solver,'IPOPT')
    params.zl = result10.zl;
    params.zu = result10.zu;
    params.lambda = result10.lambda;
end

result11 = Optimize(X0, L, U, params);