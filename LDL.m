function L=LDL(A,n)
%计算矩阵的LDL分解
%A是要分解的矩阵，n是其阶数。D的元素保存在L的对角元上，L的在相应位置。

v=zeros(n,1);
for j=1:n
    for i=1:j-1
        v(i)=A(j,i)*A(i,i);
    end
    A(j,j)=A(j,j)-A(j,1:j-1)*v(1:j-1);
    A(j+1:n,j)=(A(j+1:n,j)-A(j+1:n,1:j-1)*v(1:j-1))/A(j,j);
end
L=A;
end