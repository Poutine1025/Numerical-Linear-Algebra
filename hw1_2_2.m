%构造矩阵A和向量b
n=40;
A=hilb(n);
b=zeros(n,1);
for i=1:n
    for j=1:n
        b(i)=b(i)+1/(i+j-1);
    end
end

%用Cholesky法解方程
X=Cholesky(A,n);
L=tril(X);
x=UpperTri(L',n,LowerTri(L,n,b));

%用改进的平方根法解方程
Y=LDL(A,n);
Ly=tril(Y,-1)+eye(n);
D=triu(tril(Y));
y=UpperTri(D*Ly',n,LowerTri(Ly,n,b));

%真实解和误差
solve=linspace(1,1,40)';
error=[norm(x-solve),norm(y-solve)];