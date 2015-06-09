function result = runOPT(NSU)

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
params.NSU = 1; %number of swingups

params = getParams(params);
params.omega = zeros(params.nstates,params.N);

[X0, L, U] = getIniConBound(params);

[~, params] = conjacstructure(L, U, params);

%First result witout noise for desired trajectory
result1 = Optimize(X0, L, U, params);

% Now add more swingups and noise
params.NSU = NSU;
params = getParams(params);
params.omega = 0.001*randn(params.nstates,params.N); %Added noise, low pass filtered

[X0, L, U] = getIniConBound(params, result1);
[~, params] = conjacstructure(L, U, params);
result2 = Optimize(X0, L, U, params);

% More noise
X0 = result2.X;
params.omega = 0.01*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result3 = Optimize(X0, L, U, params);

X0 = result3.X;
params.omega = 0.05*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result4 = Optimize(X0, L, U, params);

X0 = result4.X;
params.omega = 0.1*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result5 = Optimize(X0, L, U, params);

X0 = result5.X;
params.omega = 0.5*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result6 = Optimize(X0, L, U, params);

X0 = result6.X;
params.omega = 1*randn(params.nstates,params.N); %Added noise, low pass filtered
[~, params] = conjacstructure(L, U, params);
result = Optimize(X0, L, U, params);