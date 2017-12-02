function n1=norm1(B,n)
%估计矩阵B的1范数，n为B的阶数

k=1;
x=linspace(1/n,1/n,n)';
while (k==1)
    omega=B*x; v=sign(omega); z=B'*v;
    
    if (norm(z,inf)<=z'*x)
        n1=norm(omega,1);
        k=0;
    else
        x=linspace(0,0,n)';
        j=1; z_inf=abs(z(1));
        for i=2:n
            if (abs(z(i))>z_inf)
                z_inf=abs(z(i));
                j=i;
            else
                ;
            end
        end
        x(j)=1; k=1;
    end      
end