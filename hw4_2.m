num=[0,0,0];%记录迭代次数
time=[0,0,0];%记录CPU用时
for p=1:3
    %初始化参数
    N=10*(2^p);
    h=1/N;
    %构造矩阵A和向量b
    A=zeros((N-1)^2);
    b=linspace(0,0,(N-1)^2)';
    for i=1:N-1
        for j=1:N-1
            A((i-1)*(N-1)+j,(i-1)*(N-1)+j)=4+h^2*exp(i*j*h^2);
            if j<N-1
                A((i-1)*(N-1)+j,(i-1)*(N-1)+j+1)=-1;
            end
            if j>1
                A((i-1)*(N-1)+j,(i-1)*(N-1)+j-1)=-1;
            end
        end
    end
    for k=1:(N-1)*(N-2);
        A(k,k+N-1)=-1;
        A(k+N-1,k)=-1;
    end
    for i=1:N-1
        for j=1:N-1
            b((i-1)*(N-1)+j)=h^3*(i+j);
        end
    end
    for i=1:N-1
            b((i-1)*(N-1)+1)=b((i-1)*(N-1)+1)+1;
            b((i-1)*(N-1)+N-1)=b((i-1)*(N-1)+N-1)+1;
    end
    for j=1:N-1
        b(j)=b(j)+1;
        b((N-2)*(N-1)+j)=b((N-2)*(N-1)+j)+1;
    end
    U=-triu(A,1);
    L=-tril(A,-1);
    D=A+L+U;

    %G-S迭代法
    %规定迭代停止误差
    error=0.0000001;
    M=(D-L)\U;
    g=(D-L)\b;
    y=linspace(1,1,(N-1)^2)';
    t=linspace(0,0,(N-1)^2)';
    k=0;
    T=cputime;
    while norm(y-t,inf)>=error
        t=y;
        y=M*y+g;
        k=k+1;n
    end
    time(p)=cputime-T;
    num(p)=k;
end