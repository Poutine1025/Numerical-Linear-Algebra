true_error=linspace(0,0,26);%记录真实误差
estimated_error=linspace(0,0,26);%记录估计误差
%构造矩阵A和向量b
for n=5:30
    A=eye(n);
    A(1,n)=1;
    for i=2:n
        A(i,n)=1;
        for j=1:i-1
            A(i,j)=-1;
        end
    end
    x=10*rand(n,1);
    b=A*x;
%用列主元高斯消去法解方程
    [Y,u]=GaussColumn(A,n);
    L=tril(Y,-1)+eye(n);
    U=triu(Y);
    m=b;
    for i=1:n-1
        exchange=m(i);
        m(i)=m(u(i));
        m(u(i))=exchange;
    end
    x_hat=UpperTri(U,n,LowerTri(L,n,m));
    true_error(n-4)=norm(x-x_hat,inf)/norm(x,inf);

    m=norm(b-A*x_hat,inf)/norm(b,inf);
    estimated_error(n-4)=norm1(inv(A)',n)*norm(A,inf)*m;
end
%计算真实误差的相对误差在2范数下的差值
d=norm(true_error-estimated_error);