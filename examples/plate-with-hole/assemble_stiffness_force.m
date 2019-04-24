function [KG, B1,C] = assemble_stiffness_force(num_disps, num_nodes, num_elems, ...
    cconn,knot, ngpoints1, ngpoints2, num_elem_nodes, ...
    cpts, KU, KV, m, n, k1, k2)
syms u v
E  = 1000; % aluminium
nu = 0.3;
h  = 1.0;
C = constitutive_matrix(E,nu,h);


%% Assembly
% % Assemble force
% 
% % take this as input
% tvec= zeros(2,1);
% tvec(1) = 10;
% 
% FG     = zeros(num_disps*num_nodes,1);
% for enum = 1 : num_elems
%     
%     Fe = zeros(num_disps*num_elem_nodes, 1);
%     
%     kconn     = knot.kconn;    
%     SpanRangV = knot.spanRV;          
%     elV = SpanRangV(kconn(enum,2),:); % Knot span in eta    
%     etaa = elV(1);
%     etab = elV(2);
%     [etas, wts] = gauss_quadrature_1d(ngpoints2, etaa, etab);
%         
%     % loop gauss points
%       for j = 1 : ngpoints2
%         
%         % Extract poitns and weights        
%         wt  = wts(j);
%         eta = etas(j);
%         [J1, R, ~, ~] = Jacobian_cal(cpts,KU,KV,m,n,k1,k2,cconn, enum, 0.0, eta);
%         
%         % Forming local stiffness matrix Ke
%         Fe = Fe + R*tvec*wt;
%         
%       end
%        FG(cconn(enum,:),1)= FG(cconn(enum,:),1) + Fe;
% end

% Assemble stiffness
KG     = zeros(num_disps*num_nodes, num_disps*num_nodes);
for enum = 1 : num_elems
    
    % Initialize element matrices and vectors
    Ke = zeros(num_disps*num_elem_nodes, num_disps*num_elem_nodes);
    B  = zeros(3,num_disps*num_elem_nodes);
    
    kconn     = knot.kconn;
    SpanRangU = knot.spanRU;
    SpanRangV = knot.spanRV;   
    
    elU = SpanRangU(kconn(enum,1),:); % Knot span in xi
    elV = SpanRangV(kconn(enum,2),:); % Knot span in eta
    
    % Gauss points and weights for integration
    xia  = elU(1);
    xib  = elU(2);
    etaa = elV(1);
    etab = elV(2);
    [XI, WTXI] = gauss_quadrature_2d(ngpoints1, ngpoints2, xia, xib,  etaa, etab);
    num_gauss_points = ngpoints1*ngpoints2;

    for j = 1 : num_gauss_points
        
        % Extract poitns and weights
        xivec = XI(j,:);
        wt    = WTXI(j,:);
        
        xi  = xivec(1); % Parametric coord in xi
        eta = xivec(2); % Paramtric coord in eta              
        
        %% forward coordinate transformation
        
        % xi  = 0.5*((elU(2)-elU(1)) -  xitilda + (elU(2) + elU(1))); % Eq. 44
        % eta = 0.5*((elV(2)-elV(1)) - etatilda + (elV(2) + elV(1))); % Eq. 45      
        % dxi_dxitilda   = 0.5*(elU(2) - elU(1));
        % deta_detatilda = 0.5*(elV(2) - elV(1));
        % J2 = [dxi_dxitilda, 0; 0, deta_detatilda];
        % detJ2 = det(J2);
        %% Get the first jacobian
        [J1, R1, dRdxi, dRdeta] = Jacobian_cal(cpts,KU,KV,m,n,k1,k2,cconn, enum, xi, eta);
        detJ1 = det(J1); % Eq 50
        invJ1=  inv(J1); % inverse of J1
        
        %% NURBS basis function derivatives with respect to x and y
        dRedx = [dRdxi' dRdeta']*invJ1(:,1);  %?
        dRedy = [dRdxi' dRdeta']*invJ1(:,2); %?
        
        %% Strain matrix % Equation 56
        for i2 = 1:num_elem_nodes % Loop over necp
            j1 = (2 * i2 -1);
            j2 = 2*i2;
            B(1,j1) = dRedx(i2); % 1 st row of B
            B(2,j2) = dRedy(i2); % 2 nd row of B
            B(3,j1) = dRedy(i2); % 3 rd row of B
            B(3,j2) = dRedx(i2); % 3 rd row of B
        end
        % Forming local stiffness matrix Ke
        Ke = Ke + B'*C*B*detJ1*wt;
    B1(:,:,enum)=B;
    end
    % Assemble from the local to global matricies where CnCPe corresponds
    % to the global row and column number of the global stiffness matrix
    KG(cconn(enum,:),cconn(enum,:))= KG(cconn(enum,:),cconn(enum,:)) + Ke;
end % loop elements
    
end % function