function[lamdal,lamdar,lamda]=downwardsweep(alx,blx,dlx,arx,brx,drx) %M表示初始leaf,几个参数按照upwardsweep输出输入
M=size(alx,2);
layer=binarylayer(M);
N=length(layer);
lamdal=zeros(N,M);
lamdar=lamdal;
lamda=lamdal;
lamda(N,1)=1;
for i=1:N-1
    l=layer(N-i);
    if (l==2*layer(N-i+1))
    for j=1:(l/2)
       D=2*j-1;
       E=2*j;
       B=j;
       delta=1-arx(N-i,E)*blx(N-i,D);
       lamda(N-i,D)=lamda(N-i+1,B);
       lamda(N-i,E)=lamda(N-i+1,B);
       lamdal(N-i,D)=lamdal(N-i+1,B);
       lamdar(N-i,E)=lamdar(N-i+1,B);
       lamdar(N-i,D)=lamdar(N-i+1,B)*(1-brx(N-i,E))-lamda(N-i+1,B)*drx(N-i,E)-arx(N-i,E)*lamdal(N-i+1,B)*(1-alx(N-i,D))+arx(N-i,E)*lamda(N-i+1,B)*dlx(N-i,D);
       lamdar(N-i,D)=lamdar(N-i,D)/delta;
       lamdal(N-i,E)=-blx(N-i,D)*lamdar(N-i+1,B)*(1-brx(N-i,E))+lamda(N-i+1,B)*drx(N-i,E)*blx(N-i,D)+lamdal(N-i+1,B)*(1-alx(N-i,D))-lamda(N-i+1,B)*dlx(N-i,D);
       lamdal(N-i,E)=lamdal(N-i,E)/delta;  
    end
    else
    for j=1:((l-1)/2)
       D=2*j-1;
       E=2*j;
       B=j;
       delta=1-arx(N-i,E)*blx(N-i,D);
       lamda(N-i,D)=lamda(N-i+1,B);
       lamda(N-i,E)=lamda(N-i+1,B);
       lamdal(N-i,D)=lamdal(N-i+1,B);
       lamdar(N-i,E)=lamdar(N-i+1,B);
       lamdar(N-i,D)=lamdar(N-i+1,B)*(1-brx(N-i,E))-lamda(N-i+1,B)*drx(N-i,E)-arx(N-i,E)*lamdal(N-i+1,B)*(1-alx(N-i,D))+arx(N-i,E)*lamda(N-i+1,B)*dlx(N-i,D);
       lamdar(N-i,D)=lamdar(N-i,D)/delta;
       lamdal(N-i,E)=-blx(N-i,D)*lamdar(N-i+1,B)*(1-brx(N-i,E))+lamda(N-i+1,B)*drx(N-i,E)*blx(N-i,D)+lamdal(N-i+1,B)*(1-alx(N-i,D))-lamda(N-i+1,B)*dlx(N-i,D);
       lamdal(N-i,E)=lamdal(N-i,E)/delta;  
    end
    lamda(N-i,l)=lamda(N-i+1,layer(N-i+1));
    lamdal(N-i,l)=lamdal(N-i+1,layer(N-i+1));
    lamdar(N-i,l)=lamdar(N-i+1,layer(N-i+1));
    end
end
end
