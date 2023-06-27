%初始化a=-1;c=1;M=101;sigmak=ones(2,M);timetotall=zeros(1,30);xx=zeros(1,M);yy=zeros(30,M);atl=1;
tic
[solution1,~,newsigmak,xnum]=nnnewtonit(sigmak,a,c);
xx=xnum;
yy(atl,:)=solution1;
plot(xx,solution1);
ylim([-1 1])
sigmak=newsigmak;
atl=atl+1;
toc
timetotall(atl)=toc;
