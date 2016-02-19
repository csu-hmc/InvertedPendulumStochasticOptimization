clear all
close all
clc

%Main program for pendulum swingup

% Declare parameters
params.m = 2;
params.l = 0.6;
params.g = 9.81;
params.T = 30; %Duration of motion

params.ndof  = 1;
params.nmus = 2;
params.nstates = params.ndof*2+params.nmus*2;
params.ncontrols= params.nmus;

params.NperSU = 200;
params.NSU = 1; %number of swingups
params.ineq = 0; %Enable inequality constraints

params = getParams(params);
params.muscleparam = Muscleparams();
params.solver = 'IPOPT';

[X0, L, U] = getIniConBound(params);

% Let noise increase with time
t = linspace(0,params.T,params.N);
params.omega = bsxfun(@times, 10*t/params.T, randn(params.ndof*2,params.N));
params.W = 1;

[~, params] = conjacstructure(L, U, params);

% [grad,grad_num,cjac, cjac_num] = checkDerivatives(1/2*(L+U), params); %+randn(size(X0))
% keyboard

%First result witout noise for desired trajectory
result1 = Optimize(X0, L, U, params);

[meanx,u0, uall] = plotresult(result1); 