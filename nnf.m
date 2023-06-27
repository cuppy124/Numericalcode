function[y]=nnf(x,lamda)
M=length(x);
y=zeros(M,1);
xnum=nnxnum(M);
y(1)=x(1)-exp(lamda);
y(M)=x(M)-exp(-lamda);
for i=2:M-1
   y(i)=2*((x(i+1)-x(i))/(xnum(i+1)-xnum(i))-(x(i)-x(i-1))/(xnum(i)-xnum(i-1)))/(xnum(i+1)-xnum(i-1))+2*x(i)*(x(i+1)-x(i-1))/(xnum(i+1)-xnum(i-1))*lamda*exp(lamda*xnum(i));
end
end