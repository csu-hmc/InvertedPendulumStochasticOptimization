function params = getParams(params)

params.nvarpernode = params.optstates+params.ncontrols;
params.N = params.NSU*params.NperSU;
params.h = params.T/params.NperSU;

params.nvarperSU = params.nvarpernode*params.NperSU;
params.nvars = params.nvarpernode*params.N;
params.nconSU1 = params.ncontrols*2+params.nstates*params.NperSU;
params.nconperSU = (params.ncontrols+params.optstates)*params.NperSU+params.ncontrols+params.optstates-params.nstates;
params.ncon = params.nconSU1+(params.NSU-1)*params.nconperSU+1;