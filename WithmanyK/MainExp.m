clear all
close all
clc

%Main program for pendulum swingup

% Declare parameters
params.m = 5;
params.l = 1.2;
params.g = 9.81;
params.T = 10; %Duration of motion

params.nstates  = 2;
params.ncontrols= 1;
params.optstates= 4;
params.nvarpernode = params.nstates+params.ncontrols;

params.NperSU = 60;
params.NSU = 1; %number of swingups
params = getParams(params);

[T,Y] = ode45(@(t,x) StocDynExp(t, x, params),[0 params.T],[-pi/2;0]);

plot(T,Y(:,1))
