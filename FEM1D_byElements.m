function [A,b,h,Atot]=FEM1D_byElements(x,f,alpha,kk) 
% FEM 1D Lagrangiani lineari e quadratici

% INPUT
%   x       vettore dei nodi della mesh     
%   f       handle alla funzione f(x)
%   alpha   parametro della soluzione esatta
%   kk      grado polinomiale (1=lineari, 2=quadratici)
% OUTPUT
%   A       effective system matrix
%   b       effective system rhs
%   h       vettore passi della mesh
%   Atot    system matrix

NN=length(x);       % Numero nodi
NS=NN-1;            % Numero elementi (sottointervalli)
ind=(1:NS)'; 
h=x(ind+1)-x(ind);  % h_i=x_(i+1)-x_i

switch kk 
    case 1  % LINEARI
        % Inizializzazione
        A=sparse(NN,NN); 
        b=zeros(NN,1);
        % Assembly procedure 
        for s=1:NS
            hs=h(s);                % Interval Es
            xm=(x(s)+x(s+1))/2; 
        %    
            AEl=[1 -1;-1 1]/hs;           
            % Usiamo la formula di quadratura del Punto Medio
            bEl=[1;1] * f(xm,alpha) * hs/2;
        %
            A(s:s+1,s:s+1)=A(s:s+1,s:s+1)+AEl;
            b(s:s+1)=b(s:s+1)+bEl;
        end
    case 2  % QUADRATICI
        % Inizializzazione
        A = sparse(NS*2+1,NS*2+1);
        b = zeros(NS*2+1,1);
        % Assembly procedure
        for s=1:NS
            hs=h(s);                % Interval Es
            xm=(x(s)+x(s+1))/2; 
        %    
            AEl=(1/3)*[7, -8, 1;
                       -8, 16, -8;
                       1, -8, 7];
            AEl = AEl./hs;

            x1 = x(s);
            x2 = xm;
            x3 = x(s+1);
            % Usiamo la formula di quadratura di Simpson
            bEl = (hs/6) * [f(x1, alpha); 
                              4 * f(x2, alpha); 
                              f(x3, alpha)];
        %
            idx = [2*s-1, 2*s, 2*s+1];
            A(idx,idx)=A(idx,idx)+AEl;
            b(idx)=b(idx)+bEl;
        end
end

% Extraction of the effective system matrix and rhs 
Atot=A;
A=A(2:end-1,2:end-1); 
b=b(2:end-1);
