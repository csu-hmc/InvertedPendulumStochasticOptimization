clear all
close all
clc

%Main program for pendulum swingup

% Declare parameters
params.m = 2;
params.l = 0.6;
params.g = 9.81;

params.ndof = 1;
params.nmus = 2;
params.nstates = params.ndof*2+params.nmus*2;
params.ncontrols= params.nmus;
params.muscleparam = Muscleparams();
params.solver = 'IPOPT';

params.W = 1;
% Weights = [0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100];
Times = [1 3 5 8 10 25 50 75 100 250 500 750 1000];
for i = 1:length(Times)
    for j = 1:10 %because noise
%         params.W = Weights(i);
        
        params.T = 10; %Duration of motion
        params.h = 0.2; %Time step
        params.N = params.T/params.h+1;
        params = getParams(params);
        params.omega = 0*randn(params.ndof*2,params.N); %Added noise

        [X0, L, U] = getIniConBound(params);
        [~, params] = conjacstructure(L, U, params);

        %First result witout noise for desired trajectory
        result1 = Optimize(X0, L, U, params);

        % Now add noise and increase time
        params.T = Times(i); %Duration of motion
        params.N = params.T/params.h+1;
        params = getParams(params);
        params.omega = 1*randn(params.ndof*2,params.N); %Added noise

        [X0, L, U] = getIniConBound(params);
        [~, params] = conjacstructure(L, U, params);

        result2 = Optimize(X0, L, U, params);

        filename = ['Solutions/Results_', num2str(params.T), 'time' num2str(j) '.mat'];
        save(filename, 'result1', 'result2');
    end
end