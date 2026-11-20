function [NQ, xq, wq] = quadrature_formula_1D(qf)
% INPUT
%   qf    tipo di formula di quadratura che si vuole usare
%           'simpson'   (per polinomi di grado <=3)
%           'gauss2'    (per polinomi di grado <=3)
%           'gauss3'    (per polinomi di grado <=5)
% OUTPUT
%   NQ    numero nodi della formula di quadratura
%   xq    ascisse dei nodi nel riferimento locale [0,1]
%   wq    pesi associati ai nodi per l'intervallo [0,1]

switch qf
    %---------------------------------------------------------------------          
    case 'simpson' 
        % Formula di Simpson a 3 punti su [0, 1]
        % Esatta per polinomi di grado <= 3
        NQ = 3;
        xq = zeros(1,NQ);
        wq = zeros(1,NQ);

        xq(1) = 0.0; 
        xq(2) = 0.5; % Punto medio
        xq(3) = 1.0;

        wq(1) = 1/6;  
        wq(2) = 4/6;  
        wq(3) = 1/6;

    %---------------------------------------------------------------------  
    case 'gauss2'
        % Formula di Gauss-Legendre a 2 punti su [0, 1]
        % Esatta per polinomi di grado <= 3
        NQ = 2;
        xq = zeros(1,NQ);
        wq = zeros(1,NQ);

        xq(1) = 0.5-sqrt(3)/6; 
        xq(2) = 0.5+sqrt(3)/6;

        wq(1) = 0.5;  
        wq(2) = 0.5;

    %--------------------------------------------------------------------- 
    case 'gauss3' 
        % Formula di Gauss-Legendre a 3 punti su [0, 1]
        % Esatta per grado <= 5
        NQ = 3;
        xq = zeros(1,NQ);
        wq = zeros(1,NQ);

        xq(1) = (-sqrt(3/5)+1)/2; 
        xq(2) = 0.5; 
        xq(3) = (sqrt(3/5)+1)/2;

        wq(1) = 5/18;        
        wq(2) = 8/18;        
        wq(3) = 5/18;

    %---------------------------------------------------------------------      
    otherwise    
        disp('Formula di quadratura 1D non definita!')
        return 
end

return