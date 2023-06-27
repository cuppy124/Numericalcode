function[solu]=localsolver(a,b,K) %a,b表示区间端点，K表示order,solu第一列为P^(-1)f,第二列为P^(-1)phil,第三列为P^(-1)phir,第四列为alpha(lx),beta(lx),delta(lx),第五列为alpha(rx),beta(rx),delta(rx)
xnum=zeros(1,K);
pnum=zeros(1,K);
qnum=zeros(1,K);
fnum=zeros(1,K);
glnum=zeros(1,K);
grnum=glnum;
philnum=xnum;
phirnum=xnum;
for i=1:K
 xnum(i)=(b-a)*0.5*cos((2*K-2*i+1)*pi/(2*K))+(b+a)/2;   %取chebyshev点
 pnum(i)=p(xnum(i));
 qnum(i)=q(xnum(i));
 fnum(i)=f(xnum(i));
 glnum(i)=gl(xnum(i));
 grnum(i)=gr(xnum(i));
 phirnum(i)=phir(xnum(i));
 philnum(i)=phil(xnum(i));
end
P=zeros(K,K);
for j=1:K
  P(j,j)=1;
  P(j,:)=P(j,:)+philnum(j)*(Il(a,b,K,j).*glnum)+phirnum(j)*(Ir(a,b,K,j).*grnum);
end
solu=zeros(K,6);
solu(:,1)=P\(fnum');
solu(:,2)=P\(philnum');
solu(:,3)=P\(phirnum');
solu(1,4)=chebyshevweight(a,b,K)*(glnum.*solu(:,2)')';
solu(2,4)=chebyshevweight(a,b,K)*(glnum.*solu(:,3)')';
solu(3,4)=chebyshevweight(a,b,K)*(glnum.*solu(:,1)')';
solu(1,5)=chebyshevweight(a,b,K)*(grnum.*solu(:,2)')';
solu(2,5)=chebyshevweight(a,b,K)*(grnum.*solu(:,3)')';
solu(3,5)=chebyshevweight(a,b,K)*(grnum.*solu(:,1)')';
end