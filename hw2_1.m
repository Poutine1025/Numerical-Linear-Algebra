norm_inf=linspace(0,0,16);     %��¼��ʵ�������
y=linspace(0,0,16);            %��¼���Ƶ������
for i=5:20
    H=hilb(i);
    norm_inf(i-4)=norm1(H,i);
    y(i-4)=norm(H,inf);
end
norm(y-norm_inf)               %�������