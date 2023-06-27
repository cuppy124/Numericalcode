function[y]=Jinv(x,lamda)
M=length(x);
y=zeros(M,M);
xnum=nnxnum(M);
y(1,1)=1;
y(M,M)=1;
for i=2:M-1
    y(i,i-1)=2/((xnum(i)-xnum(i-1))*(xnum(i+1)-xnum(i-1)))-2*x(i)*lamda*exp(lamda*xnum(i))/(xnum(i+1)-xnum(i-1));
    y(i,i)=-2/((xnum(i+1)-xnum(i))*(xnum(i+1)-xnum(i-1)))-2/((xnum(i)-xnum(i-1))*(xnum(i+1)-xnum(i-1)))+2*(x(i+1)-x(i-1))*lamda*exp(lamda*xnum(i))/(xnum(i+1)-xnum(i-1));
    y(i,i+1)=2/((xnum(i+1)-xnum(i))*(xnum(i+1)-xnum(i-1)))+2*x(i)*lamda*exp(lamda*xnum(i))/(xnum(i+1)-xnum(i-1));
end
y=y^(-1);
end
