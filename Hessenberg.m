function [H,Q]=Hessenberg(A)
%计算矩阵A的上Hessenberg分解
%返回修改后的矩阵和累积的正交矩阵

[m,n]=size(A);
Q=eye(n);
for k=1:n-2
    [v,beta]=house(A(k+1:n,k));
    A(k+1:n,k:n)=(eye(n-k)-beta*v*v')*A(k+1:n,k:n);
    A(1:n,k+1:n)=A(1:n,k+1:n)*(eye(n-k)-beta*v*v');
    Q=Q*blkdiag(eye(k),(eye(n-k)-beta*v*v'));
end
H=A;