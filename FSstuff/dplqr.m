function [K,S,E] = dplqr(Q, R, x, N, params)
%DPLQR Dynamic programming of LQR design for a discrete state-space system.
%   [K,S,E] =DPLQR(Q,R,x,N)  calculates the optimal gain matrices K[n]
%   such that:
%
%   For a discrete-time state-space model SYS, u[n] = -K[n]x[n] minimizes
%
%             J = 0.5*Sum {x[n]'Qx[n] + u[n]'Ru[n]}
%
%   subject to  x[n+1] = Ax[n] + Bu[n] with n = 0...N-1. Dynamics are
%   linearized
%
%   Also returned are the the solution S of the associated algebraic
%   Riccati equation and the closed-loop eigenvalues E = EIG(A-B*K).

%   written by, I. Houtzager [2006], adapted
%   Delft Center of Systems and Control

m = params.m;
l = params.l;
g = params.g;
J = m*l^2;

% backwards iteration
K = zeros(params.ncontrols,params.nstates,N);
E = zeros(params.nstates,N);
S = zeros(params.nstates,params.nstates,N);
S(:,:,N+1) = 0;
for t = N:-1:1
    theta = x(1,t);
    A = [0 1; m*g*l*sin(theta)/J 0];
    B = [0;1/J];

    sys = ss(A,B,eye(params.nstates),[0;0]);
    sys2 = c2d(sys,params.h);
    A = sys2.A;
    B = sys2.B;
    K(:,:,t) = (R + B'*S(:,:,t+1)*B)\B'*S(:,:,t+1)*A;
    S(:,:,t) = (A - B*K(:,:,t))'*S(:,:,t+1)*(A - B*K(:,:,t)) + K(:,:,t)'*R*K(:,:,t) + Q;
    if nargout >= 2
        E(:,t) = eig(A - B*K(:,:,t));
    end
end


