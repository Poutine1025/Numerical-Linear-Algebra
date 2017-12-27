function [X,P]=doubleQR(H)
%对矩阵H进行双重步位移的QR迭代
%返回修改后的矩阵和累积的正交矩阵

[~,n]=size(H);
if (n>2)
    m=n-1;
    P=eye(n);
    s=H(m,m)+H(n,n);
    t=H(m,m)*H(n,n)-H(m,n)*H(n,m);
    x=H(1,1)*H(1,1)+H(1,2)*H(2,1)-s*H(1,1)+t;
    y=H(2,1)*(H(1,1)+H(2,2)-s);
    z=H(2,1)*H(3,2);
    for k=0:n-3
        [v,beta]=house([x,y,z]');
        q=max(1,k);
        H(k+1:k+3,q:n)=(eye(3)-beta*v*v')*H(k+1:k+3,q:n);
        r=min(k+4,n);
        H(1:r,k+1:k+3)=H(1:r,k+1:k+3)*(eye(3)-beta*v*v');
        x=H(k+2,k+1);
        y=H(k+3,k+1);
        if k<(n-3)
            z=H(k+4,k+1);
        end
        P=P*blkdiag(eye(k),(eye(3)-beta*v*v'),eye(n-k-3));
    end
    [v,beta]=house([x,y]');
    H(n-1:n,n-2:n)=(eye(2)-beta*v*v')*H(n-1:n,n-2:n);
    H(1:n,n-1:n)=H(1:n,n-1:n)*(eye(2)-beta*v*v');
    X=H;
    P=P*blkdiag(eye(n-2),(eye(2)-beta*v*v')); 
else
    [P,S]=schur(H);
    X=S;
end