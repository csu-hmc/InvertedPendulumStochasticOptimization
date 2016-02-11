function result = theworhpthing(X0,L,U,params)
%-----------------------------------------------------------------------
%
% Pendulum swing-up optimization using WORHP
%
%-----------------------------------------------------------------------

% Specify number of variables and constraints
% as those values are never passed to Worhp
% it is not necessary to specify these here
% but may still be useful for initialisation
n = params.nvars;
m = params.ncon;

% Box constraints for X
XL = L;
XU = U;

% Initial estimates for X and Lambda
X = X0;
La = zeros(size(X0)); %Not sure what to put here

% Define bounds for G
GL = zeros(params.ncon,1);
GU = zeros(params.ncon,1);

% Zero initialisation for objective function f
% and constraints G
% It is important to keep G consistent during the
% later calls. size(G) = [1,m] and size(G) = [m,1]
% will lead Matlab to reallocate memory and
% therefor deny memory access for worhp
F = objfun(X0,params);
G = zeros(params.ncon,1);

% Initial estimates for Mu
Mu = zeros(params.ncon,1); %Not sure if there's a better guess

% Initialise sparsity structure of derivates
% according to the WORHP user manual (Coordinate Storage format)
DFrow = [1:1:params.nvars];
DGrow = params.DGrow;
DGcol = params.DGcol;
% HMrow = [1 2];
% HMcol = [1 2];

% Initialise structure without parameter file
% data = worhp();
% Or With parameter file
data = worhp('param.xml');

% data will contain the following fields:
% data.handle
% -> Contains adresses of C-pointers required by Worhp
% data.par
% -> Contains accessible parameters which can be altered by calling:
% --> worhp(data, 'NameOfParameter', newValue);
% data.status
% -> Contains status flag of worhp in C accessible by cnt.status
% data.scaleF
% -> Contains scaling factor for objective function, which must be used if parameter scaledObj is set to true
% data.newX
% -> Holds a flag if X has a new value. After reading the new value the flag can be reset for further tracking by calling:
% --> data.newX = false;      % Resets the flag in Matlab data structure
% --> worhp(data, data.newX); % Resets the flag in Worhp
% data.continue
% -> Holds a flag which makes Matlab continue the reverse communication loop.
% data.actions
% -> Contains fields which hold boolean values to control the reverse communication
% data.safemode
% -> Holds a flag which enables (true) or disables (false) the safety mode which can be used to test implemented examples
data.safemode = false;%true;
% After specification of dimensions and structures Worhp is called to initialise its internal data structures
% If any specific derivative structure is missing, replace the corresponding vector by []
% Obviously Using DG(HM)row = [] implies DG(HM)col = [] as well.
% e.g. worhp(data, XL, X, XU, GL, G, GU, La, Mu, [], [], [], [], []);
% Worhp will cause memory bugs if instead of X, La or Mu constant vectors are given to Worhp
% e.g. [1 2] instead of X will not work and likely break down Matlab
worhp(data, XL, X, XU, GL, G, GU, La, Mu, DFrow, DGrow, DGcol,[],[]);

% Reverse communication loop
while (data.continue)

   if (data.actions.worhp)
        worhp(data, X, La, F, G, Mu);
    end

    if (data.actions.output)
        worhp(data);
    end

    if (data.actions.f)
        F = objfun(X,params);
        worhp(data, F, 'evalF');
    end

    if (data.actions.g)
        GG = confun(X,params);
        G(:,1) = GG;
        worhp(data, G, 'evalG')
    end

    if (data.actions.df)
        DFVAL = objgrad(X,params);
        worhp(data, DFVAL, 'evalDF');
    end

    if (data.actions.dg)
        DGVAL = conjac(X,params);
        worhp(data, DGVAL, 'evalDG');
    end

    if (data.actions.fidif)
        worhp(data, X, La, F, G, Mu);
    end

end

% Call Worhp again for final output
worhp(data);

% ... or alternatively, use worhp_enum to get a meaningful enum name
worhp_enum(data.status)

% ... or check against a particular status constant
if (data.status == worhp_enum.OptimalSolution)
    disp('Optimal Solution Found (really!)')
end

% The last call of Worhp will free the memory used by Worhp
worhp(data.handle);

result.X = X;
result.params = params;
result.obj = F;
result.info = data.status;
end
