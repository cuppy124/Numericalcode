function[answ]=directmethod(n,m,a,b,A,B)  
N=n*m;
h=(b-a)/n;
xnum=zeros(N,1);    %x坐标
ynum=zeros(N,1);         %y坐标
answ=zeros(N,2);  
for i=1:n
    for j=1:m
      xnum((i-1)*m+j)=a+(i-1)*h+(1+cos((2*m-2*j+1)*pi/(2*m)))*h/2;
    end
end

H=zeros(N-1,1);    %表示间隔
for i=1:N-1
  H(i)=xnum(i+1)-xnum(i);
end

tranmatrix=zeros(N,N);      
tranf=zeros(N,1);                     
for k=1:N-2
    tranf(k)=f(xnum(k+1));
    tranmatrix(k,k+2)=2/(H(k+1)*(H(k)+H(k+1)))+p(xnum(k+1))/(H(k)+H(k+1));
    tranmatrix(k,k+1)=-2/(H(k+1)*(H(k)+H(k+1)))-2/(H(k)*(H(k)+H(k+1)))+q(xnum(k+1));
    tranmatrix(k,k)=2/(H(k)*(H(k)+H(k+1)))-p(xnum(k+1))/(H(k)+H(k+1));
end
tranmatrix(N-1,1)=A(1,1)-A(1,2)/H(1);
tranmatrix(N-1,2)=A(1,2)/H(1);
tranf(N-1)=B(1);
tranmatrix(N,N)=A(2,1)+A(2,2)/H(N-1);
tranmatrix(N,N-1)=-A(2,2)/H(N-1);
tranf(N)=B(2);

ynum=tranmatrix\tranf;

answ(:,1)=xnum;
answ(:,2)=ynum;
