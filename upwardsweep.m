function[alx,blx,dlx,arx,brx,drx]=upwardsweep(A)   %A表示6xM矩阵，行为alpha(lx),beta(lx),delta(lx),alpha(rx),beta(rx),delta(rx)
M=size(A,2);
layer=binarylayer(M);
N=length(layer);
alx=zeros(N,M);
blx=alx;
dlx=alx;
arx=alx;
brx=alx;
drx=alx;
alx(1,:)=A(1,:);
blx(1,:)=A(2,:);
dlx(1,:)=A(3,:);
arx(1,:)=A(4,:);
brx(1,:)=A(5,:);
drx(1,:)=A(6,:);
for i=2:N
l=layer(i);
if (2*l==layer(i-1))
    for j=1:l
    D=2*j-1;
    E=2*j;
    delta=1-arx(i-1,E)*blx(i-1,D);
    alx(i,j)=(1-alx(i-1,E))*(alx(i-1,D)-blx(i-1,D)*arx(i-1,E))/delta+alx(i-1,E);
    arx(i,j)=arx(i-1,E)*(1-brx(i-1,D))*(1-alx(i-1,D))/delta+arx(i-1,D);
    blx(i,j)=blx(i-1,D)*(1-brx(i-1,E))*(1-alx(i-1,E))/delta+blx(i-1,E);
    brx(i,j)=(1-brx(i-1,D))*(brx(i-1,E)-blx(i-1,D)*arx(i-1,E))/delta+brx(i-1,D);
    dlx(i,j)=(1-alx(i-1,E))*dlx(i-1,D)/delta+dlx(i-1,E)+(alx(i-1,E)-1)*blx(i-1,D)*drx(i-1,E)/delta;    %
    drx(i,j)=(1-brx(i-1,D))*drx(i-1,E)+drx(i-1,D)+(brx(i-1,D)-1)*arx(i-1,E)*dlx(i-1,D)/delta;
    end 
else
   for j=1:l-1
    D=2*j-1;
    E=2*j;
    delta=1-arx(i-1,E)*blx(i-1,D);
    alx(i,j)=(1-alx(i-1,E))*(alx(i-1,D)-blx(i-1,D)*arx(i-1,E))/delta+alx(i-1,E);
    arx(i,j)=arx(i-1,E)*(1-brx(i-1,D))*(1-alx(i-1,D))/delta+arx(i-1,D);
    blx(i,j)=blx(i-1,D)*(1-brx(i-1,E))*(1-alx(i-1,E))/delta+blx(i-1,E);
    brx(i,j)=(1-brx(i-1,D))*(brx(i-1,E)-blx(i-1,D)*arx(i-1,E))/delta+brx(i-1,D);
    dlx(i,j)=(1-alx(i-1,E))*dlx(i-1,D)/delta+dlx(i-1,E)+(alx(i-1,E)-1)*alx(i-1,D)*drx(i-1,E)/delta;
    drx(i,j)=(1-brx(i-1,D))*drx(i-1,E)+drx(i-1,D)+(brx(i-1,D)-1)*arx(i-1,E)*dlx(i-1,D)/delta;
   end
   alx(i,l)=alx(i-1,2*l-1);
   arx(i,l)=arx(i-1,2*l-1);
   blx(i,l)=blx(i-1,2*l-1);
   brx(i,l)=brx(i-1,2*l-1);
   dlx(i,l)=dlx(i-1,2*l-1);
   drx(i,l)=drx(i-1,2*l-1);
end
end
end