function [xend, rmsu, rmsu0,K,K1] = getLQRresults(result1,result2,Q,R)


params = result1.params;
params.Q = Q;
params.R = R;
X = result1.X;
params1 = result2.params;
params1.Q = Q;
params1.R = R;
X1 = result2.X;
params1.K = result2.X(end-1);
params1.Kd = result2.X(end);
params.K = result2.X(end-1);
params.Kd = result2.X(end);
if isfield(params1,'Ks')
    params.K = result2.X(end-119:end-60);
    params.Kd = result2.X(end-59:end);
    params1.K = params.K;
    params1.Kd = params.Kd;
else
    params1.Ks = 2;
    params.Ks = 2;
end


%% No feedback or noise
stdev = 0;
feedback = 0;

[x,u] = doSimulation(params, X, stdev, feedback);
[x1,u1] = doSimulation(params1, X1, stdev,feedback);

rmsu0(1) = rms(u);
rmsu0(2) = rms(u1);

%% Feedback and noise
s = rng;
stdev = 1;
feedback = 1;

[x,u,K] = doSimulation(params, X, stdev, feedback);
rng(s);
[x1,u1,K1] = doSimulation(params1, X1, stdev,feedback);

rmsu(1) = rms(u);
rmsu(2) = rms(u1);

xend(1,:) = x(:,end);
xend(2,:) = x1(:,end);
% figure
% plot(x')
% hold on
% plot(x1')