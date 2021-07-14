%=======================================================================================================
% This contain all the information for running main
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================
%
%  DATI= struct( 'name',              % set the name of the test 
%                'Domain',            % set the domain [x1,x2]
%                'bc'                 % boundary conditions DD-Dirichlet,
%                                     % NN-Neumann, PP-peridic,
%                                     % AA-absorbing
%                'c2'                 % c^2 wave speed
%                'T'                  % final time
%                'dt'                 % time step 
%                'u0',                % Initial condition u               
%                'v0',                % Initial condition  du/dt        
%                'exact_sol',         % set the exact solution
%                'force',             % set the forcing term
%                'grad_exact_1',      % set the first componenet of the gradient of the exact solution
%                'fem',               % set finite element space
%                'nqn_1D',            % number of quadrature nodes for integrals over lines
%                'refinement_vector', % set the level of refinement for the grid
%                'visual_graph',      % if you want to display the graphical results ['Y','N']
%                'print_out',         % if you want to print out the results ['Y','N']
%                'plot_errors'        % you want to print the computed errors ['Y','N']
% 
%========================================================================================================
%
%  REFRACTION
%
% Test 1 - wave propagation rightward - heterogeneous bar 
% Null initial conditions - null force
% Dirichlet conditions 
% u(0,t) = sin(2 pi t) t<= 0.5
% u(L,t) = 0 t>=0
%
%          | 1 x <= 1;
% c^2(x) = |   
%          | 2 x > 1;
%
% Example of usage:
% [errors,solutions,femregion,Dati] = C_main1D('Test1',6)
%
%*************************************************************************
%
% Test 2 - wave propagation rightward - heterogeneous bar 
% Null initial conditions - null force
% Dirichlet conditions 
% u(0,t) = sin(2 pi t) t<= 0.5
% u(L,t) = 0 t>=0
%
%          | 2 x <= 1;
% c^2(x) = |   
%          | 1 x > 1;
%
% Example of usage:
% [errors,solutions,femregion,Dati] = C_main1D('Test2',6)
%
%
%==========================================================================
%
%  GIBBS - PHENOMENON
%
% Test 6 - Triangular wave - homogeneous bar 
% Null initial conditions - null force
% Dirichlet conditions 
% u(0,t) = u(L,t) = 0 t>=0
%
% c^2(x) = 1; 
%
% Example of usage:
% [errors,solutions,femregion,Dati] = C_main1D('Test6',6)
%
%*************************************************************************
%
% Test 7 - Square wave - homogeneous bar 
% Null initial conditions - null force
% Dirichlet conditions 
% u(0,t) = u(L,t) = 0 t>=0
%
% c^2(x) = 1; 
%
% Example of usage:
% [errors,solutions,femregion,Dati] = C_main1D('Test7',6)
%
%
%==========================================================================



function [Dati]=C_dati(test)

if test=='Test1'

%dL=1/2;
c=num2str(1);
%alfa=num2str(1);
Dati = struct( 'name',             test,...
               ... % Test name
                'bc',           'DD',...                          
               ... % boundary conditions 
               'c_tilde', ['-',c,'.*(x<=0)+',c,'.*(x>=2)'],...
                ...% velocity damping layer
               'c2',             '1+0*x', ...
               ... % wave speed ...
               'kk',             '0.*x', ...
               ...  %damping coefficient
               'T',               3, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      'exp(-5.*(x-1).^2)', ...
               ... % Initial condition u               
               'v0',      '0.*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.*x.*t',...      
               ... % Definition of exact solution
               'grad_exact',     '0.*x.*t',...   
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',     '0.*t',...   
               ... % c2du/dx(0,t) 
               'neumann2',     '0.*t',...   
               ... % c2du/dx(L,t) 
               'fem',              'P3',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [4,5,6,7],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );

elseif test=='Test2' 
Dati = struct( 'name',             test,...
               ... % Test name 
               'domain',   [0 - 1/2,2 + 1/2],...                          
               ... % Domain bounds
               'sigma', '10.*(x<=0)+10.*(x>=2)',...
                ...% damping factor
                'c_tilde', '-1.*(x<=0)+1.*(x>=2)',...
                ...% velocity damping layer
               'bc',           'DD',...                          
               ... % boundary conditions                      
               'c2',             '1+0.*x', ...
               ... % wave speed ...
               'kk',             '0.*x', ...
               ...  %damping coefficient
               'T',               3, ...
               ... % Final time ...
               'dt',            0.005, ...
               ... % Time step
               'u0',      'exp(-5.*(x-1).^2)', ...
               ... % Initial condition u               
               'v0',      '0.*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.*x.*t',...      
               ... % Definition of exact solution
               'grad_exact',     '0.*x.*t',...   
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',     '0.*t',...   
               ... % c2du/dx(0,t) 
               'neumann2',     '0.*t',...   
               ... % c2du/dx(L,t) 
               'fem',              'P3',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [4,5,6,7],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );             
end