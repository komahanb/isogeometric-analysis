function [K, F] = apply_boundary_conditions(kmat, fvec, bc_dofs)
K = kmat;
F = fvec;
numbcs = size(bc_dofs);
for i = numbcs
    K(:,i) = [];
    K(:,i) = [];
    F(i) = [];
end
end