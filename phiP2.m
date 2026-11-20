function z=phiP2(i,x)
switch i
    case 1
        z = (2*x-1).*(x-1);
    case 2
        z = 4*x.*(1-x);
    case 3
        z = x.*(2*x-1);
end