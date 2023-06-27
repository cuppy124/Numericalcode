%初始化ats=1;timetotal=zeros(30,1);
tic
a=-1;
b=1;
A=zeros(2,2);
A(1,1)=1;
A(2,1)=1;
B=zeros(2,1);
B(1)=-1;
B(2)=1;
solu=directmethod(60,11,a,b,A,B);
xx=solu(:,1);
yy=solu(:,2);
plot(xx,yy)
title('60 points')
toc
timetotal(ats)=toc;
ats=ats+1;