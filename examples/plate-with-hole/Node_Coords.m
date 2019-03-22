% Used for contructing Q1 mesh
function [S] = Node_Coords(n1, p1, U1, m1, q1, V1, P, dim1, u1, v1)
uspan = findspan(n1, p1, u1, U1); % find span of knot u along xi
vspan = findspan(m1, q1, v1, V1); % find span of knot v along eta
Nu = basisfun(uspan, u1, p1, U1); % compute basis at uspan
Nv = basisfun(vspan, v1, q1, V1); % compute basis at vspan
uind = uspan - p1;                % find no. of control points along xi in e
S    = zeros(1,dim1);             % initialize the control points nodal array
for l = 0:q1                              % loop over order along eta
    temp = zeros(1,dim1);                 % Initialize a temp variable
    vind = vspan - q1 + l;                % Find no of control points along eta in e
    for i = 0:p1                          % loop over order along xi
        CP = P(uind+i+1 + vind*(n1+1),:); % get proj nodes in an element
        temp = temp + Nu(i+1)*CP;         % Update the temp for an element
    end
    S = S + Nv(l+1)*temp; % projection node array for an element
end
end