function X=Gauss(A,n)
%�ò�ѡ��Ԫ��Gauss��ȥ������A��LU�ֽ�
%AΪϵ������nΪA�Ľ��������XΪ�ֽ�

for k=1:n-1
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
    A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);  
end
X=A;
end