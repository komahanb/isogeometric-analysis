function [K, F] = apply_boundary_conditions(kmat, fvec, bc_dofs)
K = kmat;
F = fvec;
k=0;
numbcs = size(bc_dofs,2);
for i = 1: numbcs
    K(:,bc_dofs(i)-k) = [];
    K(bc_dofs(i)-k,:) = [];
    F(bc_dofs(i)-k) = [];
    k=k+1;
end
end