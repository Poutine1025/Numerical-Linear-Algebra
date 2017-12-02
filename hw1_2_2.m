%�������A������b
n=40;
A=hilb(n);
b=zeros(n,1);
for i=1:n
    for j=1:n
        b(i)=b(i)+1/(i+j-1);
    end
end

%��Cholesky���ⷽ��
X=Cholesky(A,n);
L=tril(X);
x=UpperTri(L',n,LowerTri(L,n,b));

%�øĽ���ƽ�������ⷽ��
Y=LDL(A,n);
Ly=tril(Y,-1)+eye(n);
D=triu(tril(Y));
y=UpperTri(D*Ly',n,LowerTri(Ly,n,b));

%��ʵ������
solve=linspace(1,1,40)';
error=[norm(x-solve),norm(y-solve)];