clc; clear all ;close all;

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
%surf(geometry(:,:,1),geometry(:,:,2),geometry(:,:,3))
%plot3(geometry(:,:,1),geometry(:,:,2),geometry(:,:,3))

%% 5. Connectivity Control Points
m = size(CPArray,2)-1;
n = size(CPArray,3)-1;
k1=size(KnotU,1)-1-m;
k2=size(KnotV,1)-1-n;
cconn = CPassembly_arrays(m+1,k1-1,n+1,k2-1);

%% 6. Knot Vector Connectivity array
uniqU=unique(KnotU);
uniqV=unique(KnotV);
[kconn,spanRU,spanRV]=KnotConnectivity(m+1,k1-1,n+1,k2-1,uniqU,uniqV);

% pack data into knot object
knot.kconn  = kconn;
knot.spanRU = spanRU;
knot.spanRV = spanRV;

%% 7. Jacobian of Transformation Calculation
%elem_no=2; %The element for which you want to find the jacobian
%[J, dRdxi, dRdeta]= Jacobian_cal(CPArray,KnotU,KnotV,m,n,k1,k2,cconn,elem_no)

%% 8. Assemble forces and stiffness matrix
global_conn = sort([2*cconn-1, 2*cconn]);
num_disps = 2;                  % u and v
num_nodes = (m+1)*(n+1);        % num nodes in the mesh (nodes are control points)
num_elems = (m-k1+2)*(n-k2+2);  % elements in the mesh
num_1d_gauss_points = 2;
num_elem_nodes = 9;
[FG, KG] = assemble_stiffness_force(num_disps, num_nodes, num_elems, ...
    global_conn, knot, 2, 2, num_elem_nodes, ...
    CPArray,KnotU,KnotV,m,n,k1,k2)

%% 9. Apply Boundary conditions and solve


