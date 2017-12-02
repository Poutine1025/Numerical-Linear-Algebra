%输入数据
y=[1,0.8125,0.75,1,1.3125,1.75,2.3125]';
t=[-1,-0.75,-0.5,0,0.25,0.5,0.75]';
%构造矩阵
T=ones(7,3);
for i=1:7
    T(i,1)=t(i)^2;
    T(i,2)=t(i);
end
[m,n]=size(T);
%对矩阵进行QR分解
[X,d]=QRhouse(T);
Q=eye(m);
for i=1:n
    v=[1;X(i+1:m,i)];
    Hhat=eye(m+1-i)-d(i)*v*v';
    H=blkdiag(eye(i-1),Hhat);
    Q=Q*H;
end
%求解最小二乘问题
Q1=Q(1:m,1:n);
c1=Q1'*y;
R1=Q'*T;
R=R1(1:n,1:n);
%计算的系数
coef=UpperTri(R,n,c1);