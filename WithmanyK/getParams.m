function params = getParams(params)

params.nvarpernode1 = params.nstates+params.ncontrols*3;
params.nvarpernode = params.nstates;
params.N = params.NSU*params.NperSU;
params.h = params.T/(params.NperSU-1);

params.norows = params.nvarpernode1+params.nvarpernode*(params.NSU-1);
params.nvarSU1 = params.nvarpernode1*params.NperSU;
params.nvarSU = params.nvarpernode*params.NperSU;
params.nvars = params.nvarSU1+params.nvarSU*(params.NSU-1);
params.nconSU = params.nstates*params.NperSU;
params.nconeq = params.nconSU*params.NSU+params.nstates;
params.nconineq = params.NSU;
params.ncon = params.nconeq + params.nconineq;


% if params.ineq == 1
%     params.ncon = params.ncon+params.N;
% end
