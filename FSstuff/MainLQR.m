clear all
close all
clc

load('Solutions\LinFeed\Solution_30SUs9.mat')
% load('TryForLQR.mat')
% load('TEST1.mat')
% load('WithManyK\WithManyKresult2.mat')
% result(1) = result1;
% resultuse = result7;
resultuse = result(7);

%% No feedback or noise
s = rng;

stdev = 0;
feedback = 0;

params = result(1).params;
X = result(1).X;
params.R = 1;
params.Q = 1e7;
params.K = 0;
params.Kd = 0;
if isfield(params,'Ks')
    params.K = zeros(60,1);
    params.Kd = zeros(60,1);
else
    params.Ks = 2;
end
[x,u,K,timesFS,params] = doSimulation(params, X, stdev, feedback);

% Do same thing using stochastic trajectory
rng(s);
params1 = resultuse.params;
X1 = resultuse.X;
params1.R = 1;
params1.Q = 1;
params1.K = 0;
params1.Kd = 0;
if isfield(params1,'Ks')
    params1.K = zeros(60,1);
    params1.Kd = zeros(60,1);
else
    params1.Ks = 2;
end

[x1,u1,K1,timesFS1,params1] = doSimulation(params1, X1, stdev,feedback);

figure
plot(timesFS1,x1(1,:)', 'r')
hold on
plot(timesFS1,params1.xdes(1,:)', 'r--')
plot(timesFS,x(1,:)', 'b')
plot(timesFS,params.xdes(1,:)', 'b--')
xlabel('Time [s]')
ylabel('Angle [deg]')

figure
plot(timesFS1,u1, 'r')
hold on
plot(timesFS,u, 'b')
xlabel('Time [s]')
ylabel('Control Input [Nm]')

rmsu0(1) = rms(u);
rmsu0(2) = rms(u1);

rmsx0(1) = rms(x(1,:)-params.xdes(1,:));
rmsx0(2) = rms(x1(1,:)-params1.xdes(1,:));

xend0(:,1) = x(:,end)-[pi/2;0];
xend0(:,2) = x1(:,end)-[pi/2;0];

%% Feedback no noise
s = rng;

stdev = 0;
feedback = 1;

params = result(1).params;
X = result(1).X;
params.R = 1;
params.Q = 1;
params.K = resultuse.X(end-1);
params.Kd = resultuse.X(end);
if isfield(params,'Ks')
    params.K = resultuse.X(end-119:end-60);
    params.Kd = resultuse.X(end-59:end);
else
    params.Ks = 2;
end

[x,u,K,timesFS,params] = doSimulation(params, X, stdev, feedback);

% Do same thing using stochastic trajectory
rng(s);
params1 = resultuse.params;
X1 = resultuse.X;
params1.R = 1;
params1.Q = 1;
params1.K = resultuse.X(end-1);
params1.Kd = resultuse.X(end);
if isfield(params1,'Ks')
    params1.K = resultuse.X(end-119:end-60);
    params1.Kd = resultuse.X(end-59:end);
else
    params1.Ks = 2;
end
[x1,u1,K1,timesFS1,params1] = doSimulation(params1, X1, stdev,feedback);

figure
plot(timesFS1,x1(1,:)', 'r')
hold on
plot(timesFS1,params1.xdes(1,:)', 'r--')
plot(timesFS,x(1,:)', 'b')
plot(timesFS,params.xdes(1,:)', 'b--')
xlabel('Time [s]')
ylabel('Angle [deg]')

figure
plot(timesFS1,u1, 'r')
hold on
plot(timesFS,u, 'b')
xlabel('Time [s]')
ylabel('Control Input [Nm]')

rmsunonoise(1) = rms(u);
rmsunonoise(2) = rms(u1);

rmsxnonoise(1) = rms(x(1,:)-params.xdes(1,:));
rmsxnonoise(2) = rms(x1(1,:)-params1.xdes(1,:));

xendnonoise(:,1) = x(:,end)-[pi/2;0];
xendnonoise(:,2) = x1(:,end)-[pi/2;0];


%% Feedback and noise
s = rng;

stdev = 1;
feedback = 1;

params = result(1).params;
X = result(1).X;
params.R = 1;
params.Q = 1e7;
params.K = resultuse.X(end-1);
params.Kd = resultuse.X(end);
if isfield(params,'Ks')
    params.K = resultuse.X(end-119:end-60);
    params.Kd = resultuse.X(end-59:end);
else
    params.Ks = 2;
end
[x,u,K,timesFS,params] = doSimulation(params, X, stdev, feedback);

% Do same thing using stochastic trajectory
rng(s);
params1 = resultuse.params;
X1 = resultuse.X;
params1.R = 1;
params1.Q = 1e7;
params1.K = resultuse.X(end-1);
params1.Kd = resultuse.X(end);
if isfield(params1,'Ks')
    params1.K = resultuse.X(end-119:end-60);
    params1.Kd = resultuse.X(end-59:end);
else
    params1.Ks = 2;
end
[x1,u1,K1,timesFS1,params1] = doSimulation(params1, X1, stdev,feedback);

figure
plot(timesFS1,x1(1,:)', 'r')
hold on
plot(timesFS1,params1.xdes(1,:)', 'r--')
plot(timesFS,x(1,:)', 'b')
plot(timesFS,params.xdes(1,:)', 'b--')
xlabel('Time [s]')
ylabel('Angle [deg]')

figure
plot(timesFS1,u1, 'r')
hold on
plot(timesFS,u, 'b')
xlabel('Time [s]')
ylabel('Control Input [Nm]')

rmsu(1) = rms(u);
rmsu(2) = rms(u1);

rmsx(1) = rms(x(1,:)-params.xdes(1,:));
rmsx(2) = rms(x1(1,:)-params1.xdes(1,:));

xend(:,1) = x(:,end)-[pi/2;0];
xend(:,2) = x1(:,end)-[pi/2;0];

% pendanim(x,params, 'DetSU.avi')
% pendanim(x1,params1, 'StocSU.avi')