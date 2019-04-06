function [N] = basis_function(u,knot1,n,k)
    N=zeros(n+1,1,size(u,2));
    if n==3
        for i=1:50
            N(1,1,i) = (1-2*u(i))*(1-2*u(i));
            N(2,1,i) = 2*u(i)*(1-2*u(i))+(1-u(i))*2*u(i);
            N(3,1,i) = 2*u(i)*u(i);
            N(4,1,i) = 0;
        end
        for i=51:100
            N(1,1,i) = 0;
            N(2,1,i) = 2*(1-u(i))*(1-u(i));
            N(3,1,i) = 2*u(i)*(1-u(i))+2*(1-u(i))*(2*u(i)-1);
            N(4,1,i) = (2*u(i)-1)*(2*u(i)-1);
        end
    end
    if n==2
        for i=1:100
            N(1,1,i) = (1-u(i))*(1-u(i));
            N(2,1,i) = 2*u(i)*(1-u(i));
            N(3,1,i) = u(i)*u(i);
        end
    end
end