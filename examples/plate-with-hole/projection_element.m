CPsx = geometry.coefs(1,:)'; % control points along xi direction
CPsY = geometry.coefs(2,:)'; % control points along eta direction
noCPsX = length(CPsX); % Get the number of control points along xi
noCPsY = length(CPsY); % Get the number of control points along eta

noKnotU = length(KnotU); % Number of knots in xi
noKnotV = length(KnotV); % Number of knots in eta

Cparray1 = Cparray(:,:)'; % Transpose of control points array
Dimn = size(Cparray1,2); % dimension of the control points array

Node = zeros(noKnotU*noKnotV,2); % Initial node array
count = 1; % counter variable

% building nodes for the projected Q1 mesh
for i = 1:noKnotV % loop over knots in eta
    eta = KnotV(i); % get knot i from Knot V
    for j = 1:noKnotU % Loop over Knots in xi
        xi = KnotU(j); % Get knot j from KnotU
        temp = Node_Coords(noCpsX-1, p, KnotU, noCPsY-1, q, KnotV, CParray1, Dimn, xi, eta);        
        Node(count,1) = temp(1)/temp(3); % store knot in xi
        Node(count,2) = temp(2)/temp(3); % store knot in eta
        count = count + 1; % update counter
    end
end


% Building the projection Q1 mesh
Q1_mesh = zeros(nelU*nelV,4); % Intialize the mesh array
for i = 1:nelU % loop of element in eta
    for j = 1:nelU % loop over element in xi
        Q1_mesh( (nelU*(i-1) + j), :)  = [nelU*(i-1)+j+(i-1), ...
            nelU*(i-1)+j+1+(i-1),...
            nelU*(i-1)+j+1+(i-1)+noKnotU,...
            nelU*(i-1)+j+(i-1)+noKnotU]; % Generation of Q1 mesh                           
    end
end