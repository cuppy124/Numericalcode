function[y]=nnxnum(M)
y=zeros(1,M);
for i=1:M
    y(i)=cos((2*M-2*i+1)*pi/(2*M));
end
end