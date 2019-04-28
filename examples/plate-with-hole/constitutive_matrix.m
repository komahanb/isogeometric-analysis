function [C] = constitutive_matrix(E,nu,h)
%%%%%%%INPUT%%%%%%%%%
% E - Young's modulus
% nu - Poisson's ratio
% h - thickness of the plate
%%%%%%%OUTPUT%%%%%%%%%
% C - Constitutive matrix
%%%%%%%%%%%%%%%%%%%%%%

C = zeros(3,3);
a = E/((1-nu*nu));
C(1,1) = a;
C(1,2) = a*nu;
C(2,1) = a*nu;
C(2,2) = a;
C(3,3) = 0.5*a*(1-nu);
end