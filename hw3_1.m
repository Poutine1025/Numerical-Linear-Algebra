%第一个方程
%构造矩阵A和向量b
n1=56;
A1=zeros(n1);
b1=zeros(n1,1);
for i=1:n1-1
    A1(i,i)=6;
    A1(i,i+1)=1;
    A1(i+1,i)=8;
end;
A1(n1,n1)=6;
b1(1)=7;
b1(2:n1-1)=15;
b1(n1)=14;
%QR分解
[X1,d1]=QRhouse(A1);
Q1=eye(n1);
for i=1:n1
    v=[1;X1(i+1:n1,i)];
    Hhat=eye(n1+1-i)-d1(i)*v*v';
    H=blkdiag(eye(i-1),Hhat);
    Q1=Q1*H;
end
R1=Q1'*A1;
c1=Q1'*b1;
%解得x
x1=UpperTri(R1,n1,c1);

%第二个方程
%构造矩阵A和向量b
n2=100;
A2=zeros(n2);
for i=1:n2-1
    A2(i,i)=10;
    A2(i,i+1)=1;
    A2(i+1,i)=1;
end
A2(n2,n2)=10;
b2=rand(n2,1);
%QR分解
[X2,d2]=QRhouse(A2);
Q2=eye(n2);
for i=1:n2
    v=[1;X2(i+1:n2,i)];
    Hhat=eye(n2+1-i)-d2(i)*v*v';
    H=blkdiag(eye(i-1),Hhat);
    Q2=Q2*H;
end
R2=Q2'*A2;
c2=Q2'*b2;
%解得x
x2=UpperTri(R2,n2,c2);

%第三个方程
%构造矩阵A和向量b
n3=40;
A3=hilb(n3);
b3=zeros(n3,1);
for i=1:n3
    for j=1:n3
        b3(i)=b3(i)+1/(i+j-1);
    end
end
%QR分解
[X3,d3]=QRhouse(A3);
Q3=eye(n3);
for i=1:n3
    v=[1;X3(i+1:n3,i)];
    Hhat=eye(n3+1-i)-d3(i)*v*v';
    H=blkdiag(eye(i-1),Hhat);
    Q3=Q3*H;
end
R3=Q3'*A3;
c3=Q3'*b3;
%解得x
x3=UpperTri(R3,n3,c3);