conv=zeros(1,6);
for i=1:6
    max=0;
   for j=1:101
       if (abs(yy(i,j)-yy(7,j))>max)
          max=abs(yy(i,j)-yy(7,j));
       end
   end
   conv(1,i)=log10(max);
end
plot((1:6),conv)
ylim([0,-20])