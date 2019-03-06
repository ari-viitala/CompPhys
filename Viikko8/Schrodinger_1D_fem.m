% The input x is a vector that contains the positions of the nodes
% including the endpoints where boundary conditions are enforced.
function [x,E,psi] = Schrodinger_1D_fem(x)

% Initialize arrays. In this example we use full matrices since the problem
% size is not huge. However, real computations would benefit from sparse
% matrices. Note that indecies 1 and N correspond to the nodes at the
% boundary.
N = length(x);
T = zeros(N,N);
V = zeros(N,N);
S = zeros(N,N);
H = zeros(N,N);

% Define the potential function and the two parts of the basis function at
% node i
potfun = @(s) -150*exp(-40*(s-0.25).^2) -50*exp(-10*(s-0.75).^2);
basisfun_1 = @(s,i) (s-x(i-1))/(x(i)-x(i-1)); 
basisfun_2 = @(s,i) (x(i+1)-s)/(x(i+1)-x(i));

% Loop over the basis functions and construct the matrices
for i = 2:N-1    
    %Overlap matrix S_ij of basis functions for generalized eigenvalue problem
    
    %Construct matrix for kinetic energy T
    
    %Construct matrix for potential energy V
    
end

% Construct the Hamiltonian and solve the eigenvalue problem for the
% non-zero basis functions
H = T+V;
[psi,D] = eig(H(2:N-1,2:N-1),S(2:N-1,2:N-1));

% Padd the boundary zeros for the eigenfunctions to facilitate plotting and
% extract the eigenvalues
psi = [zeros(1,N-2); psi; zeros(1,N-2)];
E = diag(D);