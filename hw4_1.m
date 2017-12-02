%初始化参数,spsilon可更改
n=100; h=1/n; a=0.5; epsilon=0.0001;
%构造矩阵A和向量b
D=-(2*epsilon+h)*eye(n-1);
H=zeros(n-1);
for i=1:n-2
    H(i,i+1)=1;
end
L=-epsilon*H'; U=-(epsilon+h)*H; A=D-L-U;
b=a*h^2*linspace(1,1,n-1)';
b(n-1)=b(n-1)-(h+epsilon);
x=h*(1:n-1)';
%计算真实解
y_true=(1-a)/(1-exp(-1/epsilon))*(1-exp(-x/epsilon))+a*x;
%规定迭代停止时的误差
error=0.0001;

%Jacobi迭代
B=D\(L+U);
g1=D\b;
y1=linspace(1,1,n-1)';
t1=linspace(0,0,n-1)';
k1=0;
while norm(y1-t1,inf)>=error
    t1=y1;
    y1=B*y1+g1;
    k1=k1+1;
end

%G-S迭代
M=(D-L)\U;
g2=(D-L)\b;
y2=linspace(1,1,n-1)';
t2=linspace(0,0,n-1)';
k2=0;
while norm(y2-t2,inf)>=error
    t2=y2;
    y2=M*y2+g2;
    k2=k2+1;
end

%SOR
B=D\(L+U);
omega=2/(1+sqrt(1-max(abs(eig(B)))));
Lomega=(D-omega*L)\((1-omega)*D+omega*U);
g3=omega*(D-omega*L)\b;
y3=linspace(1,1,n-1)';
t3=linspace(0,0,n-1)';
k3=0;
while norm(y3-t3,inf)>=error
    t3=y3;
    y3=Lomega*y3+g3;
    k3=k3+1;
end
error=[norm(y1-y_true),norm(y2-y_true),norm(y3-y_true)];
num=[k1,k2,k3];