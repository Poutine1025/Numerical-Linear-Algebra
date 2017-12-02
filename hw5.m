%初始化参数
N=64; h=1/N; k=[1,3,6,30];
D=2*eye(N-1);
U=tril(triu(ones(N-1),1),1); L=U';
M=(D-L)\U;
Y=zeros(4,5);%记录误差
%GS迭代
for i=1:4
    u=linspace(0,0,N-1)';
    for j=1:N-1
        u(j)=sin(j*pi*k(i)/N);
    end
    for m=1:5
        u=M*u;
        Y(i,m)=norm(u,inf);
    end
end
%作图
plot(Y(1,1:5));
hold on;
plot(Y(2,1:5));
hold on;
plot(Y(3,1:5));
hold on;
plot(Y(4,1:5));