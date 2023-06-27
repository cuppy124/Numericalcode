vectv=zeros(1,9);
for i=1:9
    max=0;
    for j=1:M
        tempx=abs(yy(i,j)-yy(10,j));
        if (tempx>max)
            max=tempx;
        end
        convv(1,i)=log(max);
    end
end
plot((1:9),vectv)
