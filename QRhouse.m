function [X,d]=QRhouse(A)
%��Household��������QR�ֽ�
%�����޸ĺ�ľ����beta����

[m,n]=size(A);
X=A; d=linspace(0,0,n)';
for j=1:n
    if j<m
        [v,beta]=house(X(j:m,j));
        X(j:m,j:n)=(eye(m-j+1)-beta*v*v')*X(j:m,j:n);
        d(j)=beta;
        X(j+1:m,j)=v(2:m-j+1);
    end
end