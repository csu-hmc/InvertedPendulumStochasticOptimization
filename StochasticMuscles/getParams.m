function params = getParams(params)

params.nvarpernode = params.nstates;
params.nvarSU1 = params.nvarpernode*params.N;
params.nvarSU = params.nvarSU1*params.Nsamples;
params.nvars = params.ncontrols+params.nvarSU+params.Ks;
params.nconeq = params.nstates*params.N*params.Nsamples;
params.nconineq = 0;% params.nmus*params.N*params.Nsamples;%+params.ndof;
params.ncon = params.nconeq + params.nconineq;%+1;