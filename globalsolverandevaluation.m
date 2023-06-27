%s是唯一的每个函数去不同的值
function[unum,dunum,snum]=globalsolverandevaluation(xnum,K) %xnum为1x(M+1),表示M个区间，K为order
s=2;    
M=length(xnum)-1;
layer=binarylayer(M);
N=length(layer);
pinvf=zeros(K,M);   %每列一个Bi
pinvphil=pinvf;
pinvphir=pinvf;
A=zeros(6,M);
unum=zeros(M,K);  %每行表示Bi上的u值
dunum=unum;  %表示导数
snum=zeros(1,M);  %表示所有Si
for i=1:M
   solu=localsolver(xnum(i),xnum(i+1),K);
   pinvf(:,i)=solu(:,1);
   pinvphil(:,i)=solu(:,2);
   pinvphir(:,i)=solu(:,3);
   A(1,i)=solu(1,4);
   A(2,i)=solu(2,4);
   A(3,i)=solu(3,4);
   A(4,i)=solu(1,5);
   A(5,i)=solu(2,5);
   A(6,i)=solu(3,5);
end
[alx,blx,dlx,arx,brx,drx]=upwardsweep(A);
[lamdal,lamdar,lamda]=downwardsweep(alx,blx,dlx,arx,brx,drx);
sigma=zeros(M,K);   %每行表示Bi上的取值
for i=1:M
   sigma(i,:)=pinvf(:,i)'+(lamdal(1,i)*pinvphil(:,i))'+(lamdar(1,i)*pinvphir(:,i))';
end
Jl=zeros(1,M);
Jr=zeros(1,M);
for i=1:M-1
    Jl(i+1)=Jl(i)+A(3,i)+lamdal(1,i)*A(1,i)+lamdar(1,i)*A(2,i);
    j=M+1-i;
    Jr(j-1)=Jr(j)+A(6,j)+lamdal(1,j)*A(4,j)+lamdar(1,j)*A(5,j);
end
for i=1:M
    tempx=zeros(1,K);
    tempglsigma=tempx;
    tempgrsigma=tempx;
    for j=1:K
        tempx(j)=(xnum(i+1)-xnum(i))*cos((2*K-2*j+1)*pi/(2*K))/2+(xnum(i+1)+xnum(i))/2;
        tempglsigma(j)=gl(tempx(j))*sigma(i,j);
        tempgrsigma(j)=gr(tempx(j))*sigma(i,j);
    end
    for j=1:K
        templeft=Il(xnum(i),xnum(i+1),K,j);
        tempright=Ir(xnum(i),xnum(i+1),K,j);
        aaa=Jl(i)+templeft*tempglsigma';
        bbb=Jr(i)+tempright*tempgrsigma';
        unum(i,j)=ui(tempx(j))+gr(tempx(j))/s*aaa+gl(tempx(j))/s*bbb;
        dunum(i,j)=dui(tempx(j))+dgr(tempx(j))/s*aaa+dgl(tempx(j))/s*bbb;
    end
    sigmak1=0;
    for j=1:K
        sigmak1=sigmak1+2/K*sigma(i,j)*cos((K-1)*(2*K-2*j+1)*pi/(2*K));
    end
    sigmak2=0;
    for j=1:K
        sigmak2=sigmak2+2/K*sigma(i,j)*cos((K-2)*(2*K-2*j+1)*pi/(2*K));
    end
    sigmak3=0;
    for j=1:K
        sigmak3=sigmak3+2/K*sigma(i,j)*cos((K-3)*(2*K-2*j+1)*pi/(2*K));
    end
    snum(i)=abs(sigmak2)+abs(sigmak1-sigmak3);
end
end