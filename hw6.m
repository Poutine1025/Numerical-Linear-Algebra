%输入矩阵信息
A=[2,3,4,5,6;4,4,5,6,7;0,3,6,7,8;0,0,2,8,9;0,0,0,1,10];
H=A; Q=eye(5);

num=13;%迭代次数
u=eps;%机器精度
[m,n]=size(A);
M=zeros(num,4);%记录次对角元的值
ML=zeros(num,2);%记录m和l的值
for k=1:num
   %收敛性判定
   for p=2:1:n
       if (abs(H(p,p-1))<(abs(H(p,p))+abs(H(p-1,p-1)))*u)
           H(p,p-1)=0;
       end
   end
   %找最大的m
   m=n;
   for i=n:-1:2
       if (H(i,i-1)~=0)
           if (i==2 || (i>2 && H(i-1,i-2)==0))
               Z=H(i-1:i,i-1:i);
               delta=(Z(1,1)-Z(2,2))^2+4*Z(2,1)*Z(1,2);
               if (delta>=0)
                   m=n-i;
                   break;
               else
                   
               end
           else
               m=n-i;
               break;
           end
       end
   end
   %找最小的l
   l=0;
   for j=n-m:-1:2
       if(H(j,j-1)==0)
           l=j-1;
           break;
       end
   end
   ML(k,1:2)=[m,l];
   %判断是否结束
   if (m==n)
       %输出信息
       break;
   else
       %QR迭代
       [H(l+1:n-m,l+1:n-m),P]=doubleQR(H(l+1:n-m,l+1:n-m));
       Q=Q*blkdiag(eye(l),P,eye(m));
       if l>0
           H(1:l,l+1:n-m)=H(1:l,l+1:n-m)*P;
       end
       if m>0
           H(l+1:n-m,n-m+1:n)=P'*H(l+1:n-m,n-m+1:n);
       end
   end
   M(k,1:4)=[H(2,1),H(3,2),H(4,3),H(5,4)];
end