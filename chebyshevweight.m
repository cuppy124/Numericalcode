function[y]=chebyshevweight(a,b,K)
y=(2/K)*ones(1,K);
temp=zeros(1,K);
for i=1:floor((K-1)/2)
    for j=1:K
        temp(j)=cos(i*(2*K-2*j+1)*pi/K);
    end
    temp=temp*(4/(K*(2*i-1)*(2*i+1)));
    y=y-temp;
end   
y=y*((b-a)/2);
end