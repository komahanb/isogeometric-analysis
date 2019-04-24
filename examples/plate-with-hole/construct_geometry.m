clear all ;close all;

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
colormap winter;
hold on;
plot(CPArray(1,:,1),CPArray(2,:,1),'-k*','MarkerSize',10,'Color',[1,0,0],'LineWidth',2);
hold on;
plot(CPArray(1,:,2),CPArray(2,:,2),'-k*','MarkerSize',10,'Color',[1,0,0],'LineWidth',2);
hold on;
plot(CPArray(1,:,3),CPArray(2,:,3),'-k*','MarkerSize',10,'Color',[1,0,0],'LineWidth',2);

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
[KG, B, C] = assemble_stiffness_force(num_disps, num_nodes, num_elems, ...
    global_conn, knot, 4, 4, num_elem_nodes, ...
    CPArray,KnotU,KnotV,m,n,k1,k2);
FG = zeros(24,1);
FG(17)=-100;

%% Solve applying BCs
bc_nodes = [2, 7, 10, 15, 18, 23];
[K, F] = apply_boundary_conditions(KG, FG, bc_nodes);
U = K\F;
k4=1;
for i1=1:2*(m+1)*(n+1)
    if ismember(i1,bc_nodes)
        U1(i1)=0;
    else
        U1(i1)=U(k4);
        k4=k4+1;
    end
end

%% Plotting deformed Plate
CParray2=CPArray;
k5=1;
k6=2;
for k3=1:n+1
    for k4=1:m+1
        CParray2(1,k4,k3)=CPArray(1,k4,k3)+U1(k5);
        CParray2(2,k4,k3)=CPArray(2,k4,k3)+U1(k6);
        k5=k5+2;
        k6=k6+2;
    end
end
geometry2 = nrbmak(CParray2, {KnotU, KnotV});
figure;
surf(geometry2(:,:,1),geometry2(:,:,2),geometry2(:,:,3))
hold on;
plot(CParray2(1,:,1),CParray2(2,:,1),'-k*','MarkerSize',10,'Color',[0,1,0],'LineWidth',2);
hold on;
plot(CParray2(1,:,2),CParray2(2,:,2),'-k*','MarkerSize',10,'Color',[0,1,0],'LineWidth',2);
hold on;
plot(CParray2(1,:,3),CParray2(2,:,3),'-k*','MarkerSize',10,'Color',[0,1,0],'LineWidth',2);
