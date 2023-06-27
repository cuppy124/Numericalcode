function[y]=Il(a,b,K,ii)  %ab表示端点，K表示length(ynum)，i表示到第i个chebyshev点，求左积分的左乘向量，右边为各个点的取值
y=zeros(1,K);
for i=1:K
   y(i)=(b-a)*0.5*cos((i-1)*pi*(2*K-2*ii+1)*pi/(2*K));
end
ILS=zeros(K,K);
ILS(2,1)=1;
ILS(2,3)=-0.5;
ILS(1,:)=ILS(1,:)+ILS(2,:);
for i=3:K-1
    ILS(i,i-1)=1/(2*(i-1));
    ILS(i,i+1)=-ILS(i,i-1);
    ILS(1,:)=ILS(1,:)+(-1)^(i)*ILS(i,:);
end
i=K;
ILS(i,i-1)=1/(2*(i-1));
ILS(1,:)=ILS(1,:)+(-1)^(i)*ILS(i,:);
CKF=zeros(K,K);
CKF(1,:)=0.5*ones(1,K);
for i=2:K
   for j=1:K
     CKF(i,j)=cos((i-1)*pi*(2*K-2*j+1)/(2*K));
   end
end
CKF=CKF*2/K;
y=y*ILS*CKF;
end