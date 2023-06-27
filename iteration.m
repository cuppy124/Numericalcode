function[y]=iteration(x,lamda)
y=x-Jinv(x,lamda)*nnf(x,lamda);
end