function[y]=binarylayer(M)
N=1+ceil(log2(M));
y=zeros(1,N);
y(1)=M;
for i=2:N
y(i)=ceil(y(i-1)/2);
end
end
