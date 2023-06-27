function[y]=nnG0(x,t)
Yinv1t=[1 1-t;-1 1+t];
Yinv1t=Yinv1t*0.5;
Yx=[x+1 x-1;1 1];
Jt=[-0.5 (t-1)/2;0 0];
if (t>x)
   y=Yx*Jt;
else
   y=Yx*(Yinv1t+Jt);
end
end
