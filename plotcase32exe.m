%初始化lamda=2;M=1000;x=ones(M,1);timetotal=zeros(30,1);atl=1;yy=zeros(30,M);xx=nnxnum(M);
tic
newx=iteration(x,lamda);
yy(atl,:)=newx;
x=newx;
toc
timetotal(atl,1)=toc;
atl=atl+1;
