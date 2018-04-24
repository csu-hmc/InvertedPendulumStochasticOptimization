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
params.T = 10; %Duration of motion
params.N = 600;
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
 
lengthi = 5;

Ws = [0.01 0.05 0.1 0.5 1 5 10 50 100];
stdevs = 0:0.1:1;%[0 0.001 0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.5 2 5];

for j = 1:lengthi
    for k = 1:length(Ws)
        disp(k)
        params.W = Ws(k);

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
        result(1).result(k).result(j).result = Optimize(X0, L, U, params);

        clear mex

        for i = 2:length(stdevs)
            params.warmstart = 1;

            % Now add noise and increase time
            stdev = stdevs(i)/sqrt(params.h);
            randvals = randn(params.ndof*2,params.N*params.Nsamples); %Added noise
            params.omega = stdev*randvals;

            [X0, L, U] = getIniConBound(params,0,result(i-1).result(k).result(j).result);%result1);
            [~, params] = getConjacstructure(L, U, params);

            if strcmp(params.solver,'IPOPT')
                params.zl = result(i-1).result(k).result(j).result.zl;
                params.zu = result(i-1).result(k).result(j).result.zu;
                params.lambda = result(i-1).result(k).result(j).result.lambda;
            end

            result(i).result(k).result(j).result = Optimize(X0, L, U, params);
            clear mex
        end
    end
    
    save('3103Results_better.mat','result','Ws','stdevs','lengthi')
end

