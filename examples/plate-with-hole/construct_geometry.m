%% 1. Physical Model

L = 4; % length of the plate in mm
R = 1; % radius of the hole in mm, see Fig 4

%% 2. Geometry parameters
KnotU = [0 0 0 0.5 1 1 1]'; % knot vector defining the coarsest mesh along xi direction
KnotV = [0 0 0 1 1 1]';             % knot vector defining the coarsest mesh alog eta direction
w = 0.5*(1+1/sqrt(2));             % weight of control points defining the hole

%% 3. Control point set for the coarsest mesh Eq (58)
CPArray(:,:,1) = [-R, -R*w, (1-sqrt(2))*w, 0;
                   0, (sqrt(2)-1)*w, R*w, R;
                   0, 0, 0, 0;
                   1, w, w, 1];

CPArray(:,:,2) = [-2.5, -2.5, -0.75, 0;
                   0, 0.75, 2.5, 2.5;
                   0, 0, 0, 0;
                   1, 1, 1, 1];

CPArray(:,:,3) = [-L, -L, -L, 0;
                   0, L, L, L;
                   0, 0, 0, 0;
                   1, 1, 1, 1];

%% 4. Construct the NURBS described geometry
geometry = nrbmak(CPArray, {KnotU, KnotV}); % Eq. (41)
surf(geometry(:,:,1),geometry(:,:,2),geometry(:,:,3))
%plot3(geometry(:,:,1),geometry(:,:,2),geometry(:,:,3))