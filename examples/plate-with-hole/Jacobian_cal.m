function [Jsub, R2, dRdxisub, dRdetasub] = Jacobian_cal(CP,KU,KV,m,n,k1,k2,ccp,elno,uhat,vhat)
syms u v
nel=(m-k1+2)*(n-k2+2);
    k=1;
    temp=basisfun(m+1,k1-1,KU);
    temp2=basisfun(n+1,k2-1,KV);
    M=subs(temp2,'u','v');
    N=temp(elno,:);
%     denom=sym('d');
%     denom=0;
%     for i1=1:n+1
%         for j5=1:m+1
%             denom=denom+simplify(N(j5)*M(i1)*CP(4,j5,i1));
%         end
%     end
    denom=1;
    if elno==1
        R=[N(1)*M(1)*CP(4,1,1)/denom; N(2)*M(1)*CP(4,2,1)/denom; N(3)*M(1)*CP(4,3,1)/denom;N(1)*M(2)*CP(4,1,2)/denom; N(2)*M(2)*CP(4,2,2)/denom; N(3)*M(2)*CP(4,3,2)/denom; N(1)*M(3)*CP(4,1,3)/denom; N(2)*M(3)*CP(4,2,3)/denom; N(3)*M(3)*CP(4,3,3)/denom]';
    else
        R=[N(2)*M(1)*CP(4,2,1)/denom; N(3)*M(1)*CP(4,3,1)/denom; N(4)*M(1)*CP(4,4,1)/denom;N(2)*M(2)*CP(4,2,2)/denom; N(3)*M(2)*CP(4,3,2)/denom; N(4)*M(2)*CP(4,4,2)/denom; N(2)*M(3)*CP(4,2,3)/denom; N(3)*M(3)*CP(4,3,3)/denom; N(4)*M(3)*CP(4,4,3)/denom]';
    end

%     for i=1:k2
%         for j=1:k1
%             R(k)=N(j+elno-1)*M(i);
%             k=k+1;
%         end
%     end
    dRdxi=simplify(diff(R,u));
    dRdeta=simplify(diff(R,v));
    j=elno;
    k=1;
    for i=1:k1*k2
        if j==elno+3
            j=elno;
            k=k+1;
        end
        xcp(i)=CP(1,j,k);
        ycp(i)=CP(2,j,k);
        drxcpx(i)=dRdxi(i)*xcp(i);
        drycpx(i)=dRdeta(i)*xcp(i);
        drxcpy(i)=dRdxi(i)*ycp(i);
        drycpy(i)=dRdeta(i)*ycp(i);
        j=j+1;
    end
 dxdxi=simplify(sum(drxcpx));
 dxdeta=simplify(sum(drycpx));
 dydxi=simplify(sum(drxcpy));
 dydeta=simplify(sum(drycpy));
 J(1,1)=dxdxi;
 J(1,2)=dxdeta;
 J(2,1)=dydxi;
 J(2,2)=dydeta;
 Jsub = subs(J, {u, v} , [uhat, vhat]); 
 dRdxisub = subs(dRdxi, {u, v} , [uhat, vhat]);
 dRdetasub = subs(dRdeta, {u, v} , [uhat, vhat]);
 Jsub =  double(Jsub);
 dRdxisub = double(dRdxisub);
 dRdetasub = double(dRdetasub);
 
 %% Rmatrix calculation
 R1=sym('r',[2*k1*k2,2]);
 for i1=1:2
     k=1;
     for j1=1:2*k1*k2
         if rem(j1+i1,2)~=0
            R1(j1,i1)=0;
         else
             R1(j1,i1)=R(k);
             k=k+1;
         end
     end
 end
 R2=double(subs(R1,{u,v},{uhat,vhat}));
 
end