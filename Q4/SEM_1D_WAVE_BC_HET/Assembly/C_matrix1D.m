function [M,A,C_s,D_s]=C_matrix1D(Dati,femregion)
%% [M,A] = C_matrix1D(Dati,femregion)
%==========================================================================
% Assembly of the mass matrix M and the stiffness matrix A 
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%
%    OUTPUT:
%          M           : (sparse(ndof,ndof) real) mass matrix
%          A           : (sparse(ndof,ndof) real) stiffnes matrix


addpath FESpace
addpath Assembly

fprintf('============================================================\n')
fprintf('Assembling matrices M and A ... \n');
fprintf('============================================================\n')


% connectivity infos
ndof         = femregion.ndof; % degrees of freedom
nln          = femregion.nln;  % local degrees of freedom
ne           = femregion.ne;   % number of elements
connectivity = femregion.connectivity; % connectivity matrix



% quadrature nodes and weights for integrals
[nodes_1D,w_1D] = xwlgl(nln);

% evaluation of shape bases and their derivative 
[dphiq,Grad] = basis_and_der_at_lgl(nodes_1D,nln);

% Assembly begin ...
M = sparse(ndof,ndof);  % Global Mass matrix
A = sparse(ndof,ndof);  % Global Stiffness matrix
C_s = sparse(ndof,ndof); % Global Damping matrix
D_s = sparse(ndof,ndof); % Global Damping 2 matrix

for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
      
    [BJ, pphys_1D] = C_get_Jacobian(femregion.coord(iglo,:), nodes_1D);
    % BJ        = Jacobian of the elemental map 
    % pphys_2D = vertex coordinates in the physical domain 
   
    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%
    
    % Local stiffness matrix 
    % ATT: Now c^2 has to be evaluated inside the integral!!!
    x = pphys_1D; c2 = eval(Dati.c2);
    [A_loc] = C_lap_loc(Grad,w_1D,nln,BJ,c2);

    % Assembly phase for stiffness matrix
    A(iglo,iglo) = A(iglo,iglo) + A_loc; 

    %=============================================================%
    % MASS MATRIX
    %=============================================================%
    
    % Local mass matrix 
    [M_loc] = C_mass_loc(dphiq,w_1D,nln,BJ);

    % Assembly phase for mass matrix
    M(iglo,iglo) = M(iglo,iglo) + M_loc;  
    
    %=============================================================%
    % DAMPING MATRIX C sigma
    %=============================================================%
    
    sigma = eval(Dati.sigma);
    c_tilde = eval(Dati.c_tilde);
    
    % Local damping matrix 
    [C_s_loc] = CC_lap_loc(Grad,dphiq,w_1D,nln,BJ,c_tilde,sigma);

    % Assembly phase for damping matrix
    C_s(iglo,iglo) = C_s(iglo,iglo) + C_s_loc; 
    
    %=============================================================%
    % DAMPING MATRIX D sigma
    %=============================================================%
    
    
    % Local damping matrix 
    [D_s_loc] = CD_lap_loc(dphiq,w_1D,nln,BJ,sigma);

    % Assembly phase for the damping matrix
    D_s(iglo,iglo) = D_s(iglo,iglo) + D_s_loc; 
   
   

end
