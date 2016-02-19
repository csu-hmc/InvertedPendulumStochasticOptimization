function params = getParams(params)

params.nvarpernode1 = params.nstates+params.ncontrols;
params.nvarpernode = params.nstates;
params.N = params.NSU*params.NperSU;
params.h = params.T/(params.NperSU-1);

params.nvarSU1 = params.nvarpernode1*params.NperSU;
params.nvarSU = params.nvarpernode*params.NperSU;
params.nvars = params.nvarSU1+params.nvarSU*(params.NSU-1)+2;
params.nconSU = params.nstates*params.NperSU;
params.nconeq = params.nconSU*params.NSU;%+params.ndof*2;
params.nconineq = params.nmus*params.N;
params.ncon = params.nconeq + params.nconineq;
