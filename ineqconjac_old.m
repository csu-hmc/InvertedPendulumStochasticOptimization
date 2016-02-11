function J = ineqconjac(X, params,j)

NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
nvars = params.nvars;

ix = params.nconperSU*(j-1)+(1:nstates/2);
ic = 1:nstates/2;

J = spalloc(nstates/2*NperSU,nvars,params.N);

for i = 1:NperSU
    J(ic,ix) = 1;
    ix = ix+params.nvarpernode;
    ic = ic+nstates/2;
end