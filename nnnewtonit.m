function[solution1,newphik,newsigmak,xnum]=nnnewtonit(sigmak,a,c)  %sigmak是2xM
M=size(sigmak,2);
phik=zeros(2,M);
newphik=phik;
xnum=zeros(1,M);
for i=1:M
    xnum(i)=(c-a)/2*cos((2*M-2*i+1)*pi/(2*M))+(c+a)/2;
end
weight=chebyshevweight(a,c,M);
for i=1:M
    temp=zeros(2,1);
    for j=1:M
        temp=temp+weight(j)*(nnG0(xnum(i),xnum(j))*sigmak(:,j));
    end
    phik(:,i)=temp;
end
omegak=zeros(2,2,M);
gk=zeros(2,M);
for i=1:M
    omegak(:,:,i)=-nndFphi(phik(:,i),xnum(i))-nnp0(xnum(i));
    gk(:,i)=nnp0(xnum(i))*phik(:,i)+nnF(phik(:,i),xnum(i))-sigmak(:,i);
end
P=eye(2*M);
gknew=zeros(2*M,1);
for i=1:M
    gknew(i,1)=gk(1,i);
    gknew(i+M,1)=gk(2,i);
end
for i=1:M
    temp=omegak(:,:,i);
    for j=1:M
        temp1=weight(j)*(temp*nnG0(xnum(i),xnum(j)));
        P(i,j)=P(i,j)+temp1(1,1);
        P(i,M+j)=P(i,M+j)+temp1(1,2);
        P(M+i,j)=P(M+i,j)+temp1(2,1);
        P(M+i,M+j)=P(M+i,M+j)+temp1(2,2);
    end
end
deltak=P\gknew;
newsigmak=sigmak;
for i=1:M
    newsigmak(1,i)=newsigmak(1,i)+deltak(i,1);
    newsigmak(2,i)=newsigmak(2,i)+deltak(i+M,1);
end
for i=1:M
    temp=zeros(2,1);
    for j=1:M
        temp=temp+weight(j)*(nnG0(xnum(i),xnum(j))*newsigmak(:,j));
    end
    newphik(:,i)=temp;
end
solution1=zeros(1,M);
for i=1:M
    solution1(i)=xnum(i)+newphik(1,i);
end
end