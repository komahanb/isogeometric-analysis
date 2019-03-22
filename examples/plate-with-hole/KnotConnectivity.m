% 2. Formation of knot vector connectivity array
function [KnotConnectivity, SpanRangU, SpanRangV] = KnotConnectivity(n,p,m,q)
% Inputs
% n : number of control points along xi direction
% p : order of basis function along xi direction
% m : number of control points along eta direction
% q : order of basis function along eta direction

nelU             = n-p;          % number of elements along xi
nelV             = m-q;          % number of elements along eta
nel              = nelU*nelV;    % Total number of elements

knotConnectivity = zeros(nel,2);  % knot connectivity array
SpanRangU        = zeros(nelU,2); % Knot connectvity array in xi direction
SpanRangV        = zeros(nelV,2); % Knot connectivity array in eta direction
temp = 1;
for i = 1:nelV     % loop over elements in eta
    for j = 1:nelU % loop over elements in xi
        KnotConnectivity(temp,:) = [j,i];                % Equation 62
        SpanRangU(j,:)           = [unqU(j), unqU(j+1)]; % Equation 60
        SpanRangV(i,:)           = [unqV(i), unqV(i+1)]; % Equation 61
        temp = temp + 1;
    end
end
end