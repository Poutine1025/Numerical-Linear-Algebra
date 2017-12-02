%�������A������b
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

%��ѡ��Ԫ��˹��ȥ���ⷽ��
X=Gauss(A,n);
Lx=tril(X,-1)+eye(n);
Ux=triu(X);
%xΪ�õ��Ľ�
x=UpperTri(Ux,n,LowerTri(Lx,n,b));

%����Ԫ��˹��ȥ���ⷽ��
[Y,u]=GaussColumn(A,n);
Ly=tril(Y,-1)+eye(n);
Uy=triu(Y);
m=b;
for i=1:n-1
    exchange=m(i);
    m(i)=m(u(i));
    m(u(i))=exchange;
end
%yΪ�õ��Ľ�
y=UpperTri(Uy,n,LowerTri(Ly,n,m));

norm(x-x0)
norm(y-x0)