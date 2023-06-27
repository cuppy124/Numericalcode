%初始化atspace=1;xnum=linspace(-1,1,11);xnumtotal=zeros(30,5000);timetotal=zeros(20,1);lengthtotal=zeros(30,1);DblFlag=false;
tic
K=5;
C=4;
TOL=10^(-8);
[newxnum,newDblFlag,done,unum]=meshrefinement(xnum,K,C,DblFlag,TOL);
xx=zeros(1,K*(length(xnum)-1));
yy=xx;
for i=1:(length(xnum)-1)
    for j=1:K
    xx((i-1)*K+j)=(xnum(i+1)-xnum(i))/2*cos((2*K-2*j+1)*pi/(2*K))+(xnum(i+1)+xnum(i))/2;
    yy((i-1)*K+j)=unum(i,j);
    end
end
DblFlag=newDblFlag;
plot(xx,yy)
xlim([-1 1])
toc
timetotal(atspace,1)=toc;
lengthtotal(atspace)=length(xnum);
xnumtotal(atspace,1:lengthtotal(atspace))=xnum;
xnum=newxnum;
atspace=atspace+1;
