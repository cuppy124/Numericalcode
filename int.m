function[xnum,y]=directint(a,b,n,m)
h=(b-a)/n;
N=n*m;
sigmanum=zeros(N,1);
xnum=zeros(N,1);  
for i=1:n
    for j=1:m
      xnum((i-1)*m+j)=a+(i-1)*h+(1+cos((2*m-2*j+1)*pi/(2*m)))*h/2;
    end
end
P=zeros(N,N);
bb=zeros(N,1);
weight=chebyshevweight(0,h,m);
for i=1:N
    bb(i)=fterda(xnum(i));
    P(i,i)=1;
    tempp=p(xnum(i));
    tempq=q(xnum(i));
    for i1=1:n
        for j=1:m
            tempnum=(i1-1)*m+j;
            P(i,tempnum)=P(i,tempnum)+tempp*G1(xnum(i),xnum(tempnum))*weight(j)+tempq*G0(xnum(i),xnum(tempnum))*weight(j);
        end
    end 
end
sigmanum=P\bb;
y=zeros(N,1);
for i=1:N
    y(i)=ui(xnum(i));
    for i1=1:n
        for j=1:m
            tempnum=(i1-1)*m+j;
            y(i)=weight(j)*G0(xnum(i),xnum(tempnum))*sigmanum(tempnum)+y(i);
        end
    end 
end
end
