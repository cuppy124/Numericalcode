function[y]=Ir(a,b,K,ii)  %ab表示端点，K表示length(ynum)，i表示到第i个chebyshev点，求右积分的左乘向量，右边为各个点的取值
y=zeros(1,K);
for i=1:K
   y(i)=(b-a)*0.5*cos((i-1)*pi*(2*K-2*ii+1)*pi/(2*K));
end
IRS=zeros(K,K);
IRS(2,1)=-1;
IRS(2,3)=0.5;
IRS(1,:)=IRS(1,:)-IRS(2,:);
for i=3:K-1
    IRS(i,i-1)=-1/(2*(i-1));
    IRS(i,i+1)=-IRS(i,i-1);
    IRS(1,:)=IRS(1,:)-IRS(i,:);
end
i=K;
IRS(i,i-1)=-1/(2*(i-1));
IRS(1,:)=IRS(1,:)-IRS(i,:);
CKF=zeros(K,K);
CKF(1,:)=0.5*ones(1,K);
for i=2:K
   for j=1:K
     CKF(i,j)=cos((i-1)*pi*(2*K-2*j+1)/(2*K));
   end
end
CKF=CKF*2/K;
y=y*IRS*CKF;
end