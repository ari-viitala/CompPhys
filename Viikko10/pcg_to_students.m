function [resvec, iter] = pcg_to_students(N)

% Construct here the N-by-N matrix A and the right hand side b


% It is customary to use zero as an initial guess
x = zeros(N,1);

% Construct here the preconditioner. Identity matrix corresponds to the
% case of no preconditioner
P = speye(N,N);

% Initialize vectors for the iteration. In the preconditioend version we
% need four vectors: iterate x, residual r, preconditioned residual z, 
% and search direction p
r = b - A*x;
z = P\r;
p = z;

% Set maximum iterations to 200 and the relative tolerance to 1.0e-6
maxiter = 200;
reltol = 10^(-6);
k = 0;

resvec = [];
for i=1:maxiter
    k = k + 1;
    resvec = [resvec norm(r)];
    % Compute here alpha, the next iterate x1 and the next residual r1
    
    
    % Check if the algorithm has converged
    if norm(r1)/norm(b) < reltol
        break;
    end
    % Compute here the preconditioned residual z1, beta and the next search
    % direction p1
    
    
    % Update vectors for the next iteration
    x = x1;
    r = r1;
    p = p1;
    z = z1;
end
iter = k;