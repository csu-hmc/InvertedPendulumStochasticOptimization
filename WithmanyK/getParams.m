function params = getParams(params)

params.nvarpernode = params.nstates+params.ncontrols;
params.N = params.NSU*params.NperSU;
params.h = params.T/params.NperSU;

params.nvarperSU = params.nvarpernode*params.NperSU;
params.nvars = params.nvarpernode*params.N+params.Ks;
params.nconSU1 = params.nstates*params.NperSU;
params.nconperSU = (params.ncontrols+params.nstates)*params.NperSU;
params.ncon = params.nconSU1+(params.NSU-1)*params.nconperSU+2;
