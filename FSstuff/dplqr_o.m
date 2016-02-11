function [K,S,E] = dplqr(a,b,q,r,f,n)
%DPLQR Dynamic programming of LQR design for a discrete state-space system.
%   [K,S,E] = DPLQR(SYS,Q,R,F,N) calculates the optimal gain matrices K[n]
%   such that:
%
%   For a discrete-time state-space model SYS, u[n] = -K[n]x[n] minimizes
%
%             J = 0.5*Sum {x[n]'Qx[n] + u[n]'Ru[n]} + 0.5*x[N]'Fx[N]
%
%   subject to  x[n+1] = Ax[n] + Bu[n] with n = 0...N-1.
%
%   Also returned are the the solution S of the associated algebraic
%   Riccati equation and the closed-loop eigenvalues E = EIG(A-B*K).
%
%   [K,S,E] = DPLQR(A,B,Q,R,F,N) is an equivalent syntax for LTV 
%   state-space models with dynamics  x[n+1] = A[n]x[n] + B[n]u[n].

%   written by, I. Houtzager [2006]
%   Delft Center of Systems and Control

% number of inputs
ni = nargin;

% Determine which syntax is being used
switch ni
    case 5
        if ndims(a) > 2
            error('LTI models not supported for arrays of state-space models.')
        elseif isct(a)
                error('System is not discrete. Use C2D.')
        elseif hasdelay(a),
            if isdt(a),
                % Map delay times to poles at z=0 in discrete-time case
                a = delay2z(a);
            else
                error('Not supported for continuous-time systems with delays.')
            end
        end
        [A,B,C,D,Ts] = ssdata(a);
        Q = b;
        F = r;
        R = f;
        N = r;
        array = 1;
    case 6
        A = a;
        B = b;
        Q = q;
        F = f;
        R = r;
        N = n;
        if ndims(a) > 2;
            array = 2;
        else
            array = 1;
        end
    otherwise
        error('DPLQR requires at least five or six input arguments')
end

% check dimensions and symmetry
[nax ns nna] = size(A);
[nbx nb nnb] = size(B);
[nqx nq nnq] = size(Q);
[nfx nf nnf] = size(F);
[nrx nr nnr] = size(R);
if ~isequal(nna,nnb,nnq,nnr,N);
    if array == 2
        error('The number of arrays have to be equal to N.');
    end
end
if ~isequal(nbx,nax);
    error('The A and B matrices must have the same number of rows.')
end
if ~isequal(nqx,nax,nfx,ns,nq,nf);
    error('The A, Q and F matrices must be the same size.')
end
if ~isequal(nrx,nr,nb);
    error('The R matrix must be square with as many columns as B.')
end
if ~isreal(Q) || ~isreal(R)
   error('The weight matrices Q, R must be real valued.')
end

% backwards iteration
K = zeros(nb,ns,N);
E = zeros(ns,N);
S = zeros(ns,ns,N);
S(:,:,N+1) = F;
switch array
    case 1
        for t = N:-1:1
            K(:,:,t) = (R + B'*S(:,:,t+1)*B)\B'*S(:,:,t+1)*A;
            S(:,:,t) = (A - B*K(:,:,t))'*S(:,:,t+1)*(A - B*K(:,:,t)) + K(:,:,t)'*R*K(:,:,t) + Q;
            if nargout >= 2
                E(:,t) = eig(A - B*K(:,:,t));
            end
        end
    case 2
        for t = N:-1:1
            K(:,:,t) = (R(:,:,t) + B(:,:,t)'*S(:,:,t+1)*B(:,:,t))\B(:,:,t)'*S(:,:,t+1)*A(:,:,t);
            S(:,:,t) = (A(:,:,t) - B(:,:,t)*K(:,:,t))'*S(:,:,t+1)*(A(:,:,t) - B(:,:,t)*K(:,:,t)) + K(:,:,t)'*R(:,:,t)*K(:,:,t) + Q(:,:,t);
            if nargout >= 2
                E(:,t) = eig(A(:,:,t) - B(:,:,t)*K(:,:,t));
            end
        end
end

