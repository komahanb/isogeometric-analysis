KG     = zeros(2*ncp, 2*ncp);                          % Initialize the global K matrix
FG     = zeros(2*ncp, 1);                              % Initialize the global F vector
CnnCPe = sort([2*ConnectivityCP-1, 2*ConnectivityCP]);  
noGPs  = 2*(p+1)*(q+1);                                % Total number of gauss points
[Gpwt, Gpnt] = Gauss_rule(noGPs);                      % Gauss points and weights

for i = 1:nel                                 % Loop over total number of elements
    
    Ke  = zeros(2*necp,2*necp);               % Initialize local element stiffness matrix Ke
    B   = zeros(3,2*necp);                    % Initialize strain-displacement matrix B
    elU = SpanRangU(knotConnectivity(i,1),:); % Knot span in xi
    elV = SpanRangV(knotConnectivity(i,2),:); % Knot span in eta
    
    for j = 1:noGPs % loop over total number of Gauss points
        
        gpnt     = Gpnt(j,:); % j-th Gauss point
        wg       = Gpwt(j,:); % j-th Gauss weight
        
        xitilda  = gpnt(1); % Parametric coord in xi
        etatilda = gpnt(2); % Paramtric coord in eta
        
        % coordinate transformation
        xi  = 0.5*((elU(2)-elU(1)) -  xitilda + (elU(2) + elU(1))); % Eq. 44
        eta = 0.5*((elV(2)-elV(1)) - etatilda + (elV(2) + elV(1))); % Eq. 45
        
        % jacobian of transformation of coodinates (Differentiate above
        % eqns)
        dxi_dxitilda   = 0.5*(elU(2) - elU(1));
        deta_detatilda = 0.5*(elV(2) - elV(1));
        
        % determinant of jacobian of transformation?
        J2 = dxi_dxitilda*deta_detatilda; % Eq. 49
        
        % Evaluation of NURBS basis function derivatives
        [dRedxi, dRedeta] = nrbbasisfunder({xi, eta}, geometry); % Eq 47
        
        % Jacobian matrix JI
        jacob = ...
            [geometry.coefs(1,ConnectivityCP(:,i));
            geometry.coefs(2,ConnectivityCP(:,i))]...
            *[dRedxi', dRedeta'];
        
        J1 = det(jacob); % Eq 50
        invJ1=  inv(jacob); % inverse of J1
        
        % NURBS basis function derivatives with respect to x and y
        drdx = [dRdxi'*dRdeta']*invJ1(:,1); %?
        drdy = [dRdxi'*dRdeta']*invJ1(:,2); %?
        
        % Strain matrix % Equation 56
        for i2 = 1:necp % Loop over necp
            j1 = (2 * i2 -1);
            j2 = 2*i2;
            B(1,j1) = dRedx(i2); % 1 st row of B
            B(2,j2) = dRedy(i2); % 2 nd row of B
            B(3,j1) = dRedy(i2); % 3 rd row of B
            B(3,j2) = dRedx(i2); % 3 rd row of B
        end
        
        % Forming local stiffness matrix Ke
        Ke = Ke + B'*C*B*J1*J2*Gpwt;
    end
    
    % Assemble from the local to global matricies where CnCPe corresponds
    % to the global row and column number of the global stiffness matrix
    KG(CnnCPe(i,:),CnnCPe(i,:))= KG(CnnCPe(i,:),CnnCPe(i,:)) + Ke;
end