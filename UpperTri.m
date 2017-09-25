function x=UpperTri(U,n,y)
%用回代法解上三角形方程组
%L为上三角矩阵，y为等号右边的列向量，n为U的阶数

for j=n:-1:2
    y(j)=y(j)/U(j,j);
    y(1:j-1)=y(1:j-1)-y(j)*U(1:j-1,j);
end
y(1)=y(1)/U(1,1);
x=y;
end