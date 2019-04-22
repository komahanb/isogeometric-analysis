function [Jsub, R, dRdxisub, dRdetasub] = Jacobian_cal(CP,KU,KV,m,n,k1,k2,ccp,elno, uhat, vhat)
syms u v
nel=(m-k1+2)*(n-k2+2);
%for counter=1:nel
    k=1;
    temp=basisfun(m+1,k1-1,KU);
    temp2=basisfun(n+1,k2-1,KV);
    M=subs(temp2,'u','v');
    N=temp(elno,:);
    for i=1:k2
        for j=1:k1
            R(k)=N(j)*M(i);
            k=k+1;
        end
    end
    dRdxi=simplify(diff(R,u));
    dRdeta=simplify(diff(R,v));
    j=elno;
    k=1;
    for i=1:k1*k2
        if j==m+2
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
end