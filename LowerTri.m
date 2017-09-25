function x=LowerTri(L,n,b)
%用前代法解下三角形方程组
%L为下三角矩阵，b为等号右边的列向量，n为L的阶数

for j=1:n-1
    b(j)=b(j)/L(j,j);
    b(j+1:n)=b(j+1:n)-b(j)*L(j+1:n,j);
end
b(n)=b(n)/L(n,n);
x=b;
end