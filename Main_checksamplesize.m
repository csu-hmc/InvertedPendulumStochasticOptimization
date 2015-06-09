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

% Now add more swingups and noise
NSUs = [50 100 150];% [600:100:1000];
for i = 1:length(NSUs)
    NSU = NSUs(i);
    disp(NSUs(i))
    for j = 1:10
        result = runOPT(NSU);
        filename = ['Solution_', num2str(NSUs(i)), 'SUs', num2str(j), '.mat'];
        save(filename, 'result')
    end
end
