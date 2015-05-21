function params = getParams(params)

params.nvarpernode = params.nstates+params.ncontrols;
params.N = params.NSU*params.NperSU;
params.h = params.T/params.NperSU;

params.nvarperSU = params.nvarpernode*params.NperSU;
params.nvars = params.nvarpernode*params.N+2;
params.nconSU1 = params.nvarpernode*2+params.nstates*(params.NperSU-1);
params.nconperSU = params.nvarpernode*2+(params.nstates+params.ncontrols/2)*(params.NperSU-1);
params.ncon = params.nconSU1+(params.NSU-1)*params.nconperSU;
