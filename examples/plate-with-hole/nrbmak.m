function [geometry] = nrbmak(CParray, Knots)
    KnotU = Knots{1};
    KnotV = Knots{2};
    m = size(CParray,2)-1;
    n = size(CParray,3)-1;
    k1 = size(KnotU,1)-1-m;
    k2 = size(KnotV,1)-1-n;
    u = linspace(0,KnotU(m+k1+1));
    v = linspace(0,KnotV(n+k2+1));
    N1 = zeros(m+1,1,size(u,2));
    N2 = zeros(n+1,1,size(v,2));
    %N1 = basis_function(u,KnotU,m,k1);
    %N2 = basis_function(v,KnotV,n,k2);
    for i1=1:100
        N1(:,1,i1)=basisfun(m,u(i1),k1-1,KnotU);
        N2(:,1,i1)=basisfun(n,v(i1),k2-1,KnotV);
    end
    geometry=zeros(size(u,2),size(v,2),3); %Final Points on the surface
    P1=zeros(size(u,2),size(v,2),3); %Numerator 
    P2=zeros(size(u,2),size(v,2),3); %Denominator
    for k3=1:3
        for i1 = 1:100
            for j1 = 1:100
                for i2 = 1:m+1
                    for j2 = 1:n+1
                        P1(i1,j1,k3)= P1(i1,j1,k3)+CParray(4,i2,j2)*N1(i2,1,i1)*N2(j2,1,j1)*CParray(k3,i2,j2);
                        P2(i1,j1,k3) = P2(i1,j1,k3)+CParray(4,i2,j2)*N1(i2,1,i1)*N2(j2,1,j1);
                    end
                end
                geometry(i1,j1,k3)=P1(i1,j1,k3)/P2(i1,j1,k3);
            end
        end
    end
end