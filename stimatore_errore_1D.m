function eta2 = stimatore_errore_1D(xnodes, uh, namef, alpha, kk)
% INPUT
%   xnodes      vettore nodi della mesh (estremi degli elementi)
%   uh          vettore soluzione numerica completo (inclusi i bordi = 0)
%   namef       handle alla funzione f(x)
%   alpha       parametro della soluzione esatta
%   kk          grado polinomiale (1=lineari, 2=quadratici)
% OUTPUT
%   eta2        vettore stimatore di errore locale eta^2 su ogni elemento della mesh

NS = length(xnodes)-1;      % numero di elementi (sottointervalli)
h_local = diff(xnodes);     % vettore che contiene i passi h della mesh
eta2 = zeros(NS, 1);        % inizializziamo a 0

% Inizializziamo i vettori per le derivate ai bordi di ogni elemento
uh_prime_left = zeros(NS,1);    % u'_h valutata all'inizio dell'elemento
uh_prime_right = zeros(NS,1);   % u'_h valutata alla fine dell'elemento

% CALCOLO DELLE DERIVATE DI uh SUI SINGOLI ELEMENTI
for ie = 1:NS
    hk = h_local(ie);
    switch kk
        case 1  % Caso LINEARE -> uh' è costante su tutto l'elemento
            val_derivata = (uh(ie+1)-uh(ie))/hk;
            uh_prime_left(ie) = val_derivata;
            uh_prime_right(ie) = val_derivata;

        case 2  % Caso QUADRATICO -> uh' è lineare
            % Caso Quadratico: uh'(xi) = uh_1*dphi1 + uh_2*dphi2 + uh_3*dphi3
            % I nodi locali nel riferimento [0, 1] sono x1=0, x2=0.5, x3=1
            
            ind = [2*ie-1, 2*ie, 2*ie+1];       % DoFs associati all'elemento ie
            uh_elem = uh(ind);
            
            % Valutiamo uh' all'estremo sinistro (xi=0) -> dphi=[-3, 4, -1]
            uh_prime_left(ie) = (uh_elem(1)*(-3) + uh_elem(2)*(4) + uh_elem(3)*(-1)) / hk;
            
            % Valutiamo uh' all'estremo destro (xi=1) -> dphi=[1, -4, 3]
            uh_prime_right(ie) = (uh_elem(1)*(1) + uh_elem(2)*(-4) + uh_elem(3)*(3)) / hk;
    end
end

% CALCOLO DEI SALTI (JUMP) SUI NODI INTERNI
jump = zeros(NS+1,1); 
% per le BCs jump(1)=jump(NS+1)=0
for n=2:NS
    % Calcoliamo il salto su ogni nodo interno
    jump(n) = uh_prime_left(n)-uh_prime_right(n-1);
end

% CALCOLO DI eta^2 ELEMENTO PER ELEMENTO
% Usiamo una formula di quadratura di Simpson a 3 punti per integrare il residuo interno
[~, xq_simp, wq_simp] = quadrature_formula_1D('simpson'); 

for ie = 1:NS
    hk = h_local(ie);
    x1 = xnodes(ie);
    
    % Calcolo della derivata seconda uh'' (costante sull'elemento)
    switch kk
        case 1  % Caso LINEARE
            uh_double_prime = 0;
        case 2  % Caso QUADRATICO
            ind = [2*ie-1, 2*ie, 2*ie+1];
            % phi'' nel riferimento [0, 1] è pari a [4, -8, 4]
            uh_double_prime = (uh(ind(1))*4 + uh(ind(2))*(-8) + uh(ind(3))*4) / (hk^2);
    end
    
    % Integrazione del residuo interno ||f + u''_h||^2_0,K tramite Simpson
    integrale_residuo = 0;
    for j = 1:3
        x_fisico = x1 + xq_simp(j) * hk;
        f_val = namef(x_fisico, alpha);
        
        % Termine (f(x) + u''_h(x))^2
        residuo_punto = (f_val + uh_double_prime)^2;
        integrale_residuo = integrale_residuo + wq_simp(j) * residuo_punto * hk;
    end
    
    % Termine del salto associato ai due nodi dell'elemento ie (nodo ie e nodo ie+1)
    termine_salto = 0.5 * hk * ((jump(ie))^2 + (jump(ie+1))^2);
    
    % Somma finale per l'elemento ie
    eta2(ie) = (hk^2 * integrale_residuo) + termine_salto;
end

end
