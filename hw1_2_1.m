%�������A������b
n=100;
A=zeros(n);
for i=1:n-1
    A(i,i)=10;
    A(i,i+1)=1;
    A(i+1,i)=1;
end
A(n,n)=10;
b=rand(n,1);

%��Cholesky���ⷽ��
X=Cholesky(A,n);
L=tril(X);
x=UpperTri(L',n,LowerTri(L,n,b));

%�øĽ���ƽ�������ⷽ��
Y=LDL(A,n);
Ly=tril(Y,-1)+eye(n);
D=triu(tril(Y));
y=UpperTri(D*Ly',n,LowerTri(Ly,n,b));

%��ʵ��
solve=A\b;
%���ַ��������
error=[norm(x-solve),norm(y-solve)];