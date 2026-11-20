function y=gradexactu(x,alpha)
y = alpha*x.^(alpha-1)-(alpha+1)*x.^alpha;