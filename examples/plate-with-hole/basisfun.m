function N=basisfun(n,p,U)  

syms u
unq=unique(U);
unqlen=size(unique(U),1);
B=sym('a',[n,p+1,unqlen-1]);

for i1=1:unqlen-1
    for i2=1:n+p+1
        if U(i2)==unq(i1+1)
            pos=i2-1;
            break;
        end
    end
    j1=1;
    for i2=1:pos-1 
        B(i2,j1,i1)=0;
    end
    B(pos,j1,i1)=1;
    for i2=pos+1:n
        B(i2,j1,i1)=0;
    end
    for j1=2:p+1
        for i2=1:pos-p-1
            B(i2,j1,i1)=0;
        end
        
        for i2=pos-p:pos
            if i2==0
                continue;
            end
            if  (U(i2+j1-1)-U(i2))==0 && (U(i2+j1)-U(i2+1))==0
                B(i2,j1,i1)=0;
                continue;
            
            elseif U(i2+j1-1) - U(i2) == 0
                B(i2,j1,i1)=(U(i2+j1)-u)/(U(i2+j1)-U(i2+1))*B(i2+1,j1-1,i1);
                continue
                
            elseif U(i2+j1)- U(i2+1) == 0
                B(i2,j1,i1)=(u-U(i2))/(U(i2+j1-1)-U(i2))*B(i2,j1-1,i1);
                continue;
                
            else
            B(i2,j1,i1)=(u-U(i2))/(U(i2+j1-1)-U(i2))*B(i2,j1-1,i1)+(U(i2+j1)-u)/(U(i2+j1)-U(i2+1))*B(i2+1,j1-1,i1);
            end
            
        end
        
        for i2=pos+1:n
            B(i2,j1,i1)=0;
        end
    end
    
end
    N=sym('n',[unqlen-1,n]);
    for i2=1:unqlen-1
        N(i2,:)=B(:,p+1,i2);
    end
end

