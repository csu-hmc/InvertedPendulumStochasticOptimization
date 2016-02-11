function params = getParams(params)

params.nvarpernode = params.nstates+params.ncontrols;
params.N = params.NSU*params.NperSU;
params.h = params.T/(params.NperSU-1);

params.nvarperSU = params.nvarpernode*params.NperSU;
params.nvars = params.nvarpernode*params.N+2;
params.nconSU1 = params.nstates*params.NperSU+params.nstates;
params.nconperSU = (params.ncontrols+params.nstates)*params.NperSU+params.nstates;
params.ncon = params.nconSU1+(params.NSU-1)*params.nconperSU+2;

if params.ineq == 1
    params.ncon = params.ncon+params.N;
end

% load('Destraj.mat');
% x = reshape(X(1:end-2),params.nvarpernode, params.NperSU);
% params.Xdes = x(1,:);