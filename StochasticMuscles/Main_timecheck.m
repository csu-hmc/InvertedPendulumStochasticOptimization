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
params.ncontrols= 0;%params.nmus;

params.muscleparam = Muscleparams();
params.solver = 'IPOPT';

params.h = 0.01; %Use value just found
params.W = 0.1;

Ts = [1 5 10 20 30 40 50 100 200];
obj = zeros(size(Ts));

for i = 1:length(Ts)
    for j = 1:5
        params.T = Ts(i);
        params.N = params.T/params.h+1;
        %Get a integer value
        params.N = round(params.N);
        params.h = params.T/(params.N-1);
        
        params = getParams(params);

        params.omega = 0*randn(params.ndof*2,params.N); %Added noise

        [X0, L, U] = getIniConBound(params);
        [~, params] = conjacstructure(L, U, params);

        %First result witout noise for desired trajectory
        result1 = Optimize(X0, L, U, params);

        % Now add noise and increase time
        params = getParams(params);
        params.omega = 60*sqrt(params.h)*randn(params.ndof*2,params.N); %Added noise

        [X0, L, U] = getIniConBound(params);
        [~, params] = conjacstructure(L, U, params);

        result2 = Optimize(X0, L, U, params);

        filename = ['Solutions/Results_', num2str(params.T), 'timecheck' num2str(j) '.mat'];
        save(filename, 'result1', 'result2');
        objs(i,j) = result2.obj;
    end
end

plot(objs)