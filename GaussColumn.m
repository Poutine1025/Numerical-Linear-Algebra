function [X,u]=GaussColumn(A,n)
%������ԪGauss��ȥ������A��LU�ֽ�
%AΪϵ������nΪA�Ľ��������XΪ�ֽ�

L=0;           %��¼��Ԫ���ڵ�����
M=0;           %��¼����Ԫ
u=zeros(n-1,1);%��¼�û�����
for k=1:n-1
    L=k;
    M=A(k,k);
    for j=k+1:n
        if (abs(A(j,k))>abs(M))
            M=A(j,k);
            L=j;
        else
            ;
        end
    end
    exchange=A(k,1:n);
    A(k,1:n)=A(L,1:n);
    A(L,1:n)=exchange;
    u(k)=L;
    if (A(k,k)~=0)
        A(k+1:n,k)=A(k+1:n,k)/A(k,k);
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
    else
        break;
    end
end
X=A;
end
