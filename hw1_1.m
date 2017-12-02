%构造矩阵A和向量b
n=84;
A=zeros(n);
b=zeros(n,1);
for i=1:n-1
    A(i,i)=6;
    A(i,i+1)=1;
    A(i+1,i)=8;
end;
A(n,n)=6;
b(1)=7;
b(2:n-1)=15;
b(n)=14;
x0=linspace(1,1,n)';

%不选主元高斯消去法解方程
X=Gauss(A,n);
Lx=tril(X,-1)+eye(n);
Ux=triu(X);
%x为得到的解
x=UpperTri(Ux,n,LowerTri(Lx,n,b));

%列主元高斯消去法解方程
[Y,u]=GaussColumn(A,n);
Ly=tril(Y,-1)+eye(n);
Uy=triu(Y);
m=b;
for i=1:n-1
    exchange=m(i);
    m(i)=m(u(i));
    m(u(i))=exchange;
end
%y为得到的解
y=UpperTri(Uy,n,LowerTri(Ly,n,m));

norm(x-x0)
norm(y-x0)