
function [Phi, Qd] = compute_system(omega, Q, dt)
A = [ -vcross(omega) vcross(omega)
       zeros(3,3)    zeros(3,3)   ];
I3 = eye(3);
% H = [I3, I3];                       % measurement
% Ob = obsv(A,H);
% unob = length(A)-rank(Ob)
[Phi] = compute_discrete(A, dt);
B = eye(6);
Qd = B*Q*B'*dt + 1/2*A*B*Q*B'*(dt)^2 + 1/2*B*Q*B'*A'*(dt)^2 + 1/3*A*B*Q*B'*A'*(dt)^3;

% compute discrete time system
function [Phi] = compute_discrete(A, dt)
Phi = eye(size(A)) + A*dt + (A*dt)^2/2 + (A*dt)^3/6;
