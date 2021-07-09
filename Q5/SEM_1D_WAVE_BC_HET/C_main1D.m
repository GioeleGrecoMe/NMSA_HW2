function [errors,solutions,femregion,Dati] = C_main1D(TestName,nRef)
%==========================================================================
% Solution of the Wave Equation with linear finite elements
% (non homogeneous Dirichlet boundary conditions)
%==========================================================================
%
%    INPUT:
%          TestName    : (string)  see C_dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Dati        : (struct)  see C_dati.m
%
% Usage:
%    [errors,solutions,femregion,Dati] = C_main1D('Test1',3)



addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing
addpath SemLib


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Dati = C_dati(TestName);
Dati.nRefinement = nRef;

%==========================================================================
% MESH GENERATION
%==========================================================================

%[Region] = C_create_mesh(Dati);
[Region] = C_create_mesh_sem(Dati);


%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = C_create_femregion(Dati,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES
%==========================================================================

[M_nbc,A_nbc] = C_matrix1D(Dati,femregion);



%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================
Dati.t = 0;
[b_nbc] = C_rhs1D(Dati,femregion);

t = 0;
if (strcmp(Dati.bc,'NN'))
    b_nbc(1)   = b_nbc(1)   - eval(Dati.neumann1);
    b_nbc(end) = b_nbc(end) + eval(Dati.neumann2);
end
%==========================================================================
% BUILD INITIAL CONDITIONS
%==========================================================================
x = femregion.coord;
u0 = eval(Dati.u0);
v0 = eval(Dati.v0);

%% First step of leapfrog

% Solve the reduced system
b_nbc = (M_nbc-0.5*Dati.dt^2*A_nbc)*u0 + Dati.dt*M_nbc*v0  + 0.5*Dati.dt^2*b_nbc;

%==========================================================================
% COMPUTE BOUNDARY CONDITIONS -- MODIFICATION OF M an b
%==========================================================================

if (strcmp(Dati.bc,'NN'))
    u1 = M_nbc\b_nbc;
elseif(strcmp(Dati.bc,'DD'))
    [M,b,u_g] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati);
    u1 = M\b;
    u1 = u1 + u_g;
elseif(strcmp(Dati.bc,'AA'))
    b_nbc(1)   = b_nbc(1) -  sqrt(eval(Dati.c2))*v0(1)*Dati.dt^2;
    b_nbc(end) = b_nbc(end) - sqrt(eval(Dati.c2))*v0(end)*Dati.dt^2;
    u1 = M_nbc\b_nbc;
else
    M_nbc(1,:) = M_nbc(1,:) + M_nbc(end,:);
    b_nbc(1)   = b_nbc(1) + b_nbc(end);
    M_nbc(end,:) = zeros(1,size(M_nbc,1));
    M_nbc(end,1) = 1;
    M_nbc(end,end) = -1;
    b_nbc(end) = 0;
    
    u1 = M_nbc\b_nbc;
    
end

[u1] = C_snapshot_1D(femregion,u1,Dati);


fprintf('============================================================\n')
fprintf('Starting time-loop ... \n');
fprintf('============================================================\n')
time_surf = [0: Dati.dt : Dati.T]; 
u_surf = zeros(length(time_surf),size(u1,1));
u_surf(1,:) = u0;
u_surf(2,:) = u1;

k_surf = 3;
for t = Dati.dt : Dati.dt : Dati.T - Dati.dt
    
    fprintf('time = %5.3e \n',t);
    
    Dati.t = t;
    [b_nbc] = C_rhs1D(Dati,femregion);
    if (strcmp(Dati.bc,'NN'))
        b_nbc(1)   = b_nbc(1)   - eval(Dati.neumann1);
        b_nbc(end) = b_nbc(end) + eval(Dati.neumann2);
    end
    
    
    b_nbc = (2*M_nbc - Dati.dt^2*A_nbc)*u1 - M_nbc*u0  + Dati.dt^2 * b_nbc;

    if (strcmp(Dati.bc,'NN'))
        u2 = M_nbc\b_nbc;
    elseif(strcmp(Dati.bc,'PP'))
        M_nbc(1,:) = M_nbc(1,:) + M_nbc(end,:);
        b_nbc(1)   = b_nbc(1) + b_nbc(end);
        M_nbc(end,:) = zeros(1,size(M_nbc,1));
        M_nbc(end,1) = 1;
        M_nbc(end,end) = -1;
        b_nbc(end) = 0;      
        u2 = M_nbc\b_nbc;
    elseif(strcmp(Dati.bc,'AA'))
        b_nbc(1)   = b_nbc(1) -  sqrt(eval(Dati.c2))*(u1(1)-u0(1))*Dati.dt;
        b_nbc(end) = b_nbc(end) - sqrt(eval(Dati.c2))*(u1(end)-u0(end))*Dati.dt;
        u2 = M_nbc\b_nbc;
    else
        [M,b,u_g] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati);
        u2 = M\b;
        u2 = u2 + u_g;
    end
   [u2] = C_snapshot_1D(femregion,u2,Dati);
    
   u_surf(k_surf,:) = u2;
   k_surf = k_surf + 1;
    % update
    u0 = u1;
    u1 = u2;
    
end
figure(100);
surf(u_surf,'EdgeColor','None');
gf = gca;
xlim([1  length(u1)]); ylim([1 Dati.T/Dati.dt]);
gf.XTick = [1 (length(u1)+1)/2 length(u1)];
gf.XTickLabel = {num2str(min(x)),num2str(round(mean(x))) ,num2str(max(x))};
gf.YTick = [0 Dati.T/2/Dati.dt Dati.T/Dati.dt];
gf.YTickLabel = {num2str(0),num2str(Dati.T*0.5) ,num2str(Dati.T)};
view(2); xlabel('space-axis'); ylabel('time-axis'); title('u_h(x,t)');
%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[solutions] = C_postprocessing(Dati,femregion,u2);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
errors = [];
if (Dati.plot_errors)
    [errors] = C_compute_errors(Dati,femregion,solutions);
end



