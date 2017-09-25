function L=Cholesky(A,n)
%计算对称正定矩阵的Cholesky分解
%A是所要分解的正定矩阵，n是其阶数。L元素的元素保存在相应位置。

for k=1:n
    A(k,k)=sqrt(A(k,k));
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
    for j=k+1:n
        A(j:n,j)=A(j:n,j)-A(j:n,k)*A(j,k);
    end
end
L=A;
end