clear all
close all
clc

%Main program for pendulum swingup

% Declare parameters
params.m = 2;
params.l = 0.6;
params.g = 9.81;

params.ndof  = 1;
params.nmus = 2;
params.nstates = params.ndof*2+params.nmus*2;
params.ncontrols= params.nmus;

params.muscleparam = Muscleparams();
params.solver = 'IPOPT';

params.T = 10; %Duration of motion
params.W = 0.1;

Ns = [100 200 300 400 500 600 700 800 900 1000];
obj = zeros(size(Ns));

for i = 1:length(Ns)
    for j = 1:5
        params.N = Ns(i);
        params.h = params.T/(params.N-1);
        params = getParams(params);

        params.omega = 0*randn(params.ndof*2,params.N); %Added noise

        [X0, L, U] = getIniConBound(params);
        [~, params] = conjacstructure(L, U, params);

        %First result witout noise for desired trajectory
        result1 = Optimize(X0, L, U, params);

        % Now add noise and increase time
        params = getParams(params);
        params.omega = 2/sqrt(params.h)*randn(params.ndof*2,params.N); %Added noise

        [X0, L, U] = getIniConBound(params);
        [~, params] = conjacstructure(L, U, params);

        result2 = Optimize(X0, L, U, params);

        filename = ['Solutions/Results_', num2str(params.N), 'hcheck10_' num2str(j) '.mat'];
        save(filename, 'result1', 'result2');
        objs(i,j) = result2.obj;
    end
end

plot(mean(objs,2))