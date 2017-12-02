function [v,beta]=house(x)
%计算Householder变换
%返回v和beta值

n=length(x);
v=linspace(0,0,n)';
eta=norm(x,inf); x=x/eta;
sigma=x(2:n)'*x(2:n);
v(2:n)=x(2:n);
if sigma==0
    beta=0;
else
    alpha=sqrt(x(1)^2+sigma);
    if x(1)<=0
        v(1)=x(1)-alpha;
    else
        v(1)=-sigma/(x(1)+alpha);
    end
    beta=2*v(1)^2/(sigma+v(1)^2); v=v/v(1);
end