function [D_loc]=CD_lap_loc(dphiq,w_1D,nln,BJ,sigma)
%% [M_loc]=C_d_loc(dphiq,w_1D,nln,BJ)
%==========================================================================
% Build the local damping matrix for the term (uv)
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          dphiq       : (array real) evaluation of the basis function on
%                        quadrature nodes
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : Jacobian of the map 
%
%    OUTPUT:
%          M_loc       :  (array real) Local mass matrix

D_loc=zeros(nln,nln);

for i=1:nln
    for j=1:nln
        for k=1:length(w_1D)
            Binv = 1./BJ;      % inverse
            Jdet = BJ;      % determinant 
            D_loc(i,j) = D_loc(i,j) + (Jdet.*w_1D(k).* sigma(k)) .* dphiq(1,k,i).* dphiq(1,k,j); %UGUALE A MASS_LOC MA MOLTIPLICATO PER SIGMA
        end
    end
end                                         
                                              

