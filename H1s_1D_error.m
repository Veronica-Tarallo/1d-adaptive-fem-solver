function H1s_error = H1s_1D_error(xnodes, gradphi, ne, qf, graduexact, uh, alpha)
% INPUT
%   xnodes      vettore nodi della mesh (estremi degli elementi)
%   gradphi     handle della funzione phi'(x)
%   ne          numero di DoFs per elemento
%   qf          tipo di formula di quadratura
%   graduexact  handle alla derivata esatta u'(x)
%   uh          vettore soluzione numerica completo (inclusi i bordi = 0)
%   alpha       parametro della soluzione esatta
% OUTPUT
%   H1s_error   errore in seminorma H1

% Estraiamo i punti e i pesi di quadratura sul riferimento 1D [0,1]
[NQ, xq, wq] = quadrature_formula_1D(qf);  

% Calcoliamo i gradienti delle funzioni di base locali nei punti di quadratura
for k = 1:ne
    ind_phi = ['ind', num2str(k)];
    gphiQ.(ind_phi) = gradphi(k, xq);
end

squared_H1s = 0;
NS = length(xnodes)-1;  % Numero di elementi (sottointervalli)

for ie = 1:NS
  % Selezione dei Gradi di Libertà globali associati all'elemento ie
  if ne == 2        % [Caso LINEARE (kk=1)]
      ind = [ie, ie+1];     % nodi dell'elemento ie sono ie, ie+1
  elseif ne == 3    % [Caso QUADRATICO (kk=2)]
      ind = [2*ie-1, 2*ie, 2*ie+1];     % nodi dell'elemento ie sono 2*ie-1, 2*ie, 2*ie+1
  end
  
  x1 = xnodes(ie);
  x2 = xnodes(ie+1);
  hs = x2-x1;       % lunghezza dell'elemento fisico
  
  % Mappatura lineare dal riferimento [0, 1] all'elemento fisico [x1, x2]
  % x(xi) = x1 + xi * hs
  xqK = x1 + xq * hs;  % nodi di quadratura sull'elemento fisico
  
  % Derivata della soluzione esatta nei punti di quadratura dell'elemento ie
  graduexact_hat = graduexact(xqK, alpha);
  
  % Derivata della soluzione discreta nel riferimento locale
  graduh_hat = zeros(1, NQ);
  for s = 1:ne 
      ind_phi_s = ['ind', num2str(s)];
      graduh_hat = graduh_hat + uh(ind(s)) * gphiQ.(ind_phi_s);  
  end 
  
  % Trasformazione della derivata dall'elemento di riferimento all'elemento fisico
  graduh_fisico = graduh_hat / hs;

  % Calcolo della differenza tra derivata esatta e numerica
  du = abs(graduexact_hat - graduh_fisico);
  squared_H1s = squared_H1s + hs*sum(wq.*(du.^2)); 
end

H1s_error = sqrt(squared_H1s);