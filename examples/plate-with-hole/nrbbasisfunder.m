function [dredxi,dredeta]=nrbbasisfunder(xi,eta)
    dredxi(1)=-4*(1-2*xi)*(1-eta)^2;
    dredxi(2)=0;
    dredxi(3)=(4-12*xi)*(1-eta)^2;
    dredxi(4)=0;
    dredxi(5)=(4*xi)*(1-eta)^2;
    dredxi(6)=0;
    dredxi(7)=-4*(1-2*xi)*(2*eta-2*eta^2);
    dredxi(8)=0;
    dredxi(9)=(4-12*xi)*(2*eta-2*eta^2);
    dredxi(10)=0;
    dredxi(11)=(4*xi)*(2*eta-2*eta^2);
    dredxi(12)=0;
    dredxi(13)=-4*(1-2*xi)*(eta^2);
    dredxi(14)=0;
    dredxi(15)=(4-12*xi)*(eta^2);
    dredxi(16)=0;
    dredxi(17)=(4*xi)*(eta^2);
    dredxi(18)=0;
    
    
    dredeta(1)=-2*(1-eta)*(1-2*xi)^2;
    dredeta(2)=0;
    dredeta(3)=(4*xi-6*xi^2)*-2*(1-eta);
    dredeta(4)=0;
    dredeta(5)=(2*xi^2)*-2*(1-eta);
    dredeta(6)=0;
    dredeta(7)=(1-2*xi)^2*(2-4*eta);
    dredeta(8)=0;
    dredeta(9)=(4*xi-6*xi^2)*(2-4*eta);
    dredeta(10)=0;
    dredeta(11)=(2*xi^2)*(2-4*eta);
    dredeta(12)=0;
    dredeta(13)=(1-2*xi)^2*(2*eta);
    dredeta(14)=0;
    dredeta(15)=(4*xi-6*xi^2)*(2*eta);
    dredeta(16)=0;
    dredeta(17)=(2*xi^2)*(2*eta);
    dredeta(18)=0;
end