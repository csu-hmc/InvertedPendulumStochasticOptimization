clear all
close all
clc

%Main program for pendulum swingup

% Declare parameters
params.m = 5;
params.l = 1.2;
params.g = 9.81;
params.ktheta = 0.1;
params.kthetadot = 1;
params.T = 10; %Duration of motion
params.K = -2; %Feedback gain 

params.nstates  = 4;
params.ncontrols= 2;
params.nvarpernode = params.nstates+params.ncontrols;

params.NperSU = 60;
params.NSU = 5; %number of swingups
params.N = params.NSU*params.NperSU;
params.h = params.T/params.NperSU;

params.nvarperSU = params.nvarpernode*params.NperSU;
params.nvars = params.nvarpernode*params.N;
params.nconSU1 = params.nvarpernode*2+params.nstates*(params.NperSU-1);
params.nconperSU = params.nvarpernode*2+(params.nstates+params.ncontrols/2)*(params.NperSU-1);
params.ncon = params.nconSU1+(params.NSU-1)*params.nconperSU;

omega = randn(params.nstates,params.N); %Added noise, low pass filtered
% a = 1;
% b = 1/20*ones(20,1);
[b,a] = butter(2,0.1);
params.omega = filter(b,a,omega);
% plot(params.omega')

[T,Y] = ode45(@(t,x) StocDynExp(t, x, params),[0 10],[0;0;0;0]);

plot(T,Y(:,1))
hold on;
plot(T,Y(:,3), 'r')