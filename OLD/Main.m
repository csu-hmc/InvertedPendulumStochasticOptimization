clear all
close all
clc

%Main program for pendulum swingup

% Declare parameters
params.m = 5;
params.l = 1.2;
params.g = 9.81;
params.ktheta = 1;
params.kthetadot = 1;
params.T = 10; %Duration of motion

params.nstates  =4;
params.ncontrols=1;
params.nvarpernode = params.nstates+params.ncontrols;

params.NperSU = 60;
params.NSU = 2; %number of swingups
params.N = params.NSU*params.NperSU;
params.h = params.T/params.NperSU;

params.nvarperSU = params.nvarpernode*params.NperSU;
params.nvars = params.nvarpernode*params.N;
params.ncon = 794;%params.nstates*params.NperSU+params.NSU*(2*params.NperSU+2*(params.NperSU-1));

params.Omega_y = 10; %Expected variance of omega
params.Omega_m = 1; %Not sure about this
params.omega = 1*randn(params.nstates,params.N); %Added noise

[X0, L, U] = getIniConBound(params);

[~, params] = conjacstructure(L, U, params);
[grad,grad_num,cjac, cjac_num] = checkDerivatives(L, U, params);

% result = Optimize(X0, L, U, params);
% plotresult(result, params);