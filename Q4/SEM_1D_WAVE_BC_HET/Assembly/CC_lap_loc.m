function [C_loc] = CC_lap_loc(Grad,dphiq,w_1D,nln,BJ,c_tilde,sigma)
%% [K_loc] = C_lap_loc(Grad,w_1D,nln,BJ)
%==========================================================================
% Build the local stiffness matrix for the term grad(u)grad(v)
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          Grad        : (array real) evaluation of the gradient on
%                        quadrature nodes
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : (array real) Jacobian of the map 
%          c2          : (array real) evaluation of wavespeed
%
%    OUTPUT:
%          K_loc       :  (array real) Local stiffness matrix


C_loc = zeros(nln,nln);

%% General implementation -- to be used with general finite element spaces
for i=1:nln
    for j=1:nln
        for k=1:length(w_1D)
            Binv = 1./BJ;    % inverse
            Jdet = BJ;       % determinant 
            C_loc(i,j) = C_loc(i,j) + (Jdet.*w_1D(k).*c_tilde(k).*sigma(k)) .* ((Grad(k,:,i) * Binv) .* dphiq(1,k,j));
        end
    end
end
%STIFFNESS MATRIX CON invece che grad*grad, grad*dphi




                                              
                                              

