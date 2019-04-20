function [geometry] = nrbmak(CParray, Knots)
    syms u
    KnotU = Knots{1};
    KnotV = Knots{2};
    unq1=unique(KnotU);
    unq2=unique(KnotV);
    unqlen1=size(unq1,1);
    unqlen2=size(unq2,1);
    m = size(CParray,2)-1;
    n = size(CParray,3)-1;
    k1 = size(KnotU,1)-1-m;
    k2 = size(KnotV,1)-1-n;
    u1 = linspace(0,KnotU(m+k1+1));
    v1 = linspace(0,KnotV(n+k2+1));
    temp1 = sym('n1',[unqlen1-1,m+1]);
    temp2 = sym('n2',[unqlen2-1,n+1]);
    temp1 = basisfun(m+1,k1-1,KnotU);
    temp2 = basisfun(n+1,k2-1,KnotV);
    N1 = zeros(m+1,size(u1,2));
    N2 = zeros(n+1,size(v1,2));
    j1=2;
    j2=2;
    for i1=1:size(u1,2)
        
        if u1(i1)<unq1(j1)
            N1(:,i1)=subs(temp1(j1-1,:),u,u1(i1));
        elseif u1(i1)==1
            N1(m+1,i1)=1;
        else
            N1(:,i1)=subs(temp1(j1,:),u,u1(i1));
            j1=j1+1;
        end
        if v1(i1)<unq2(j2)
            N2(:,i1)=subs(temp2(j2-1,:),u,v1(i1));
        elseif v1(i1)==1
            N2(n+1,i1)=1;
        else
            N2(:,i1)=subs(temp2(j2,:),u,v1(i1));
            j2=j2+1;
        end
    end
    geometry=zeros(size(u1,2),size(v1,2),3); %Final Points on the surface
    P1=zeros(size(u1,2),size(v1,2),3); %Numerator 
    P2=zeros(size(u1,2),size(v1,2),3); %Denominator
    for k3=1:3
        for i1 = 1:100
            for j1 = 1:100
                for i2 = 1:m+1
                    for j2 = 1:n+1
                        P1(i1,j1,k3)= P1(i1,j1,k3)+CParray(4,i2,j2)*N1(i2,i1)*N2(j2,j1)*CParray(k3,i2,j2);
                        P2(i1,j1,k3) = P2(i1,j1,k3)+CParray(4,i2,j2)*N1(i2,i1)*N2(j2,j1);
                    end
                end
                geometry(i1,j1,k3)=P1(i1,j1,k3)/P2(i1,j1,k3);
            end
        end
    end
    
end