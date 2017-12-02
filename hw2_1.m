norm_inf=linspace(0,0,16);     %记录真实的无穷范数
y=linspace(0,0,16);            %记录估计的无穷范数
for i=5:20
    H=hilb(i);
    norm_inf(i-4)=norm1(H,i);
    y(i-4)=norm(H,inf);
end
norm(y-norm_inf)               %计算误差