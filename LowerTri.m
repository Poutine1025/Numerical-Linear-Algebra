function x=LowerTri(L,n,b)
%��ǰ�������������η�����
%LΪ�����Ǿ���bΪ�Ⱥ��ұߵ���������nΪL�Ľ���

for j=1:n-1
    b(j)=b(j)/L(j,j);
    b(j+1:n)=b(j+1:n)-b(j)*L(j+1:n,j);
end
b(n)=b(n)/L(n,n);
x=b;
end