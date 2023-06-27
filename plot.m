length1=lengthtotal(1:10,1);
total1=sum(length1);
at1=0;
xx2=zeros(1,total1);
yy2=xx2;
for i=1:10
    for j=1:length1(i)
    xx2(at1+j)=xnumtotal(i,j);
    yy2(at1+j)=i;
    end
    at1=at1+length1(i);
end
scatter(xx2,yy2,15)
