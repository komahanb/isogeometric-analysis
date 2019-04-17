function [N] = basis_function(i,knot1,n)
    syms u v
    %N=zeros(n+1,1);
    if n==3
        if i==1
            N(1,1) = (1-2*u)*(1-2*u);
            N(2,1) = 2*u*(1-2*u)+(1-u)*2*u;
            N(3,1) = 2*u*u;
            N(4,1) = 0;
        else
            N(1,1) = 0;
            N(2,1) = 2*(1-u)*(1-u);
            N(3,1) = 2*u*(1-u)+2*(1-u)*(2*u-1);
            N(4,1) = (2*u-1)*(2*u-1);
        end
    end
    if n==2
            N(1,1) = (1-v)*(1-v);
            N(2,1) = 2*v*(1-v);
            N(3,1) = v*v;
    end
end