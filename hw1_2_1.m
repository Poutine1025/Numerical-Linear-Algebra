%构造矩阵A和向量b
n=100;
A=zeros(n);
for i=1:n-1
    A(i,i)=10;
    A(i,i+1)=1;
    A(i+1,i)=1;
end
A(n,n)=10;
b=rand(n,1);

%用Cholesky法解方程
X=Cholesky(A,n);
L=tril(X);
x=UpperTri(L',n,LowerTri(L,n,b));

%用改进的平方根法解方程
Y=LDL(A,n);
Ly=tril(Y,-1)+eye(n);
D=triu(tril(Y));
y=UpperTri(D*Ly',n,LowerTri(Ly,n,b));

%真实解
solve=A\b;
%两种方法的误差
error=[norm(x-solve),norm(y-solve)];