%Non-linear iteration and eigenvalue solving is provided.

%Three things need to be implemented for (a)
%(1) Calculation of the effective potential V_eff = VHa + Vext + Vxc.
%(2) Assemble Hamiltonian matrix H from kinetic part and effective potential.
%(3) Calculate density from eigenfunctions of H. Density is calculated from 
% 	 *normalized* eigenfunctions, which are occupied (number of occupied
%	 eigenfunctions is given by Nocc).

%For (b) you will need to add calculation of total energy.

function [psi, dens, ene] = DFT_solver(N,Nocc,nuclei)

% Set up the grid
L = pi;
h = L/N;
x = h*[1:N-1]; %Note that this discretization does *not* contain the boundary points, which are zero.

% Kinetic energy matrix T. Use your preferred method of discretization. Remember to use overlap matrix S
% appropriately if you are using non-orthogonal basis functions
T = zeros(N-1, N-1);

% Precalculate external potential Vext from the nuclei.
% (this does not change when we find self-consistent solution).
% nuclei(1) contains the position of the first nuclei etc.
Vext = zeros(size(x));


% Parameters for the non-linear iteration
mix = 0.5;
maxiter = 100;
tol = 0.001;

% Build initial H
% To get the initial guess we can construct hamiltonian via
% "H = T + Vext"
H = zeros(N-1,N-1);

% Solve initial eigenvalues and -vectors. Remember to use S in case of non-orthogonal basis functions.
[X,D] = eig(H);
[evs, perm] = sort(diag(D));
psi = X(:,perm); %Sort eigenvectors according to the eigenvalues.

% Normalize eigenfunctions and compute initial density
dens_init = zeros(N-1,1);
%Compute the density from Nocc and eigenfunctions X of Hamiltonian H

VHa = zeros(N-1,1); %electrostatic potential
dens = dens_init;

% Non-linear iteration starts here
for i=1:maxiter
    % Compute electrostatic potential VHa from the current density

    % Compute exchange-correlation potential Vxc from the current density
    Vxc = zeros(N-1,1);
    % Construct effective potential
    Veff = Vext + VHa + Vxc;
    % Build current H  and solve the eigenproblem
    H = zeros(N-1, N-1);
    
    [X,D] = eig(H);
    [evs, perm] = sort(diag(D));
    psi = X(:,perm);
    % Update density
    dens_old = dens;
    
    %Compute new density
    dens = zeros(N-1,1);
    
    % Mix old and new density. Without this kind of underrelaxation the
    % non-linear iteration won't converge
    dens = (1-mix)*dens_old + mix*dens;
    % Check for convergence
    norm(dens - dens_old)
    if sum(abs(dens - dens_old)) < tol
        break
    end
end

% For (b) you need compute final total energy
% in order to find energy minimum
ene = 0.0
