%初始化a=-1;b=1;n=20;m=5;ats=1;timetotal=zeros(30,1);
tic
[xx,yy]=directint(a,b,n,m);
plot(xx,yy)
title('100 points')
toc
timetotal(ats)=toc;
ats=ats+1;