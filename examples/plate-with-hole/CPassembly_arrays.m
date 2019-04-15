% 1. Formation of the control point assembly array
function [ConnectivityCP] = CPassembly_arrays(n,p,m,q)
% Inputs
% n : number of control points along xi direction
% p : order of basis function along xi direction
% m : number of control points along eta direction
% q : order of basis function along eta direction

nel            = (n-p)*(m-q);     % Total number of elements in the mesh
ncp            = n*m;             % Total number of control points
necp           = (p+1)*(q+1);     % Total number of control points in an element
ConnectivityCP = zeros(necp, nel); % control point assembly array

A = 0;                % loop counter
e = 0;                % loop counter
for j = 1:m           % loop over CPs along eta
    for i = 1:n       % loop over CPs along xi
        A = A+1;
        if i >= (p+1) && j >= (q+1)
            e = e +1;
            for jloc = 0:q                     % loop over order eta
                for iloc = 0:p                 % loop over order xi
                    B = A -jloc*n -iloc;      % assigning global no.
                    b = jloc*(p+1) + iloc +1;  % assigning local no.
                    ConnectivityCP(b,e) = B;   % connectivity array
                end
            end
        end
    end
end
ConnectivityCP = sort(ConnectivityCP)';
end
