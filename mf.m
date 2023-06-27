function[newxnum,newDblFlag,done,unum]=meshrefinement(xnum,K,C,DblFlag,TOL)    %done表示是否结束所有循环
[unum,~,snum]=globalsolverandevaluation(xnum,K);
up=0;
done=false;
down=0;
M=length(xnum)-1;
for i=1:M-1
    temp1=abs(unum(i+1,1)-unum(i,1));
    temp2=abs(unum(i+1,1)+unum(i,1));
    if (temp1>up)
        up=temp1;
    end
    if (temp2>down)
        down=temp2;
    end
end
TEST=up/down;
if (TEST>TOL)
change=zeros(1,M);
sdiv=max(snum)/2^(C);
for i=1:M
    if (snum(i)>=sdiv)
        change(i)=2;
    elseif (mod(2,i)==1)
        if ((snum(i)+snum(i+1)*2^(K)<sdiv))
            change(i)=-1;
            change(i+1)=-1;
        end
    end
end
newxnum=zeros(1,M+sum(change)/2);
atwhere=1;
for i=1:M
if (change(i)==0)
    newxnum(atwhere)=xnum(i);
    newxnum(atwhere+1)=xnum(i+1);
    atwhere=atwhere+1;
elseif (change(i)==2)
    newxnum(atwhere)=xnum(i);
    newxnum(atwhere+1)=(xnum(i)+xnum(i+1))/2;
    newxnum(atwhere+2)=xnum(i+1);
    atwhere=atwhere+2;
elseif (mod(2,i)==1)
    newxnum(atwhere)=xnum(i);
    newxnum(atwhere+1)=xnum(i+2);
    atwhere=atwhere+1;
end
end
newDblFlag=false;
elseif (DblFlag==false)
newxnum=zeros(1,2*M+1);
newxnum(2*M+1)=xnum(M+1);
for j=1:M
    newxnum(2*j-1)=xnum(j);
    newxnum(2*j)=(xnum(j)+xnum(j+1))/2;
end
newDblFlag=true;
else
newDblFlag=true;
done=true;
end
end
