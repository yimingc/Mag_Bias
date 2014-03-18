
function [Phi, Qd] = compute_system(omega, Q, dt)
A = [ -vcross(omega) vcross(omega)
       zeros(3,3)    zeros(3,3)   ];
[Phi] = compute_discrete(A, dt);
B = eye(6);
Qd = B*Q*B'*dt;

% compute discrete time system
function [Phi] = compute_discrete(A, dt)
Phi = eye(size(A)) + A*dt + (A*dt)^2/2 + (A*dt)^3/6;
