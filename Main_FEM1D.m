%% -----    PDE: PROGETTO 3     ----- %%
clear; clc; close all;

nameexactu=@exactu;
namegradexactu=@gradexactu;
namef=@effe;
alpha = 4/3;        % valore di alpha DA MODIFICARE [alpha = 4/3, 10]
kk = 1;             % grado dei polinomi DA MODIFICARE [k=1,2]

xmin=0; xmax=1;

switch kk
    case 1  % LINEARI
        ne = 2;
        %phi=@phiP1;
        gradphi=@gradphiP1;
        quadrature_type = 'gauss2';
    case 2  % QUADRATICI
        ne = 3;
        %phi=@phiP2;
        gradphi=@gradphiP2;
        quadrature_type = 'gauss3';
    otherwise
        disp('Error')
        return
end   

%% MESH UNIFORMI

% Creiamo la successione delle mesh
n_mesh = 9;                             % numero di mesh nella successione
nele_initalMesh = 4;                    % numero di elementi nella mesh iniziale (4 elementi)
nnodes_vector = zeros(n_mesh,1);        % inizializziamo il vettore del numero di nodi
errH1s_vettore = zeros(n_mesh, 1);      % inizializziamo il vettore degli errori

% numero di nodi per le mesh nella successione
for i=1:n_mesh
    nnodes_vector(i) = nele_initalMesh*2^(i-1)+1;
end

for i=1:n_mesh
  nodes_number = nnodes_vector(i);
  h_mesh=(xmax-xmin)/(nodes_number-1);
  xnodes=(xmin:h_mesh:xmax)';

  % Vettore che contiene tutti i DoFs
  switch kk
      case 1
          all_nodes = xnodes;
      case 2
          all_nodes = zeros((nodes_number-1)*2+1, 1);
          for j=1:length(all_nodes)
              if mod(j,2)==1
                  % inserisco i veri nodi (xnodes)
                  all_nodes(j) = xnodes((j+1)/2);
              else
                  % inserisco i nodi aggiuntivi (i punti medi)
                  all_nodes(j) = (xnodes(j/2)+xnodes(j/2+1))/2;
              end
          end
  end

  [A,b,h_vettore,Atot]=FEM1D_byElements(xnodes,namef,alpha,kk);

%
% Solution by direct method PA=LU
  N=size(A,1); 
  ua=A\b;
  ue=nameexactu(all_nodes(2:end-1),alpha); %xmin and xmax not considered

% ERRORE IN SEMINORMA 1
ua_completo = [0; ua; 0];   % ricostruiamo il vettore ai bordi
errH1s_vettore(i) = H1s_1D_error(xnodes, gradphi, ne, quadrature_type, namegradexactu, ua_completo, alpha);
end

%% METODO ADATTATIVO

errH1s_adapt_vettore = zeros(n_mesh, 1);    % inizializziamo il vettore degli errori
dof_adapt = zeros(n_mesh, 1);               % inizializziamo il vettore dei DoFs

h_initial_mesh = (xmax-xmin)/nele_initalMesh;   % passo h della mesh iniziale uniforme di 4 elementi
xnodes_adapt=(xmin:h_initial_mesh:xmax)';       % mesh iniziale (uniforme di 4 elementi)

for i = 1:n_mesh
    NS = length(xnodes_adapt) - 1;      % numero di elementi della mesh (sottointervalli)
    h_local = diff(xnodes_adapt);       % vettore che contiene gli h della mesh
    
    % 1. SOLVE
    [A,b] = FEM1D_byElements(xnodes_adapt,namef,alpha,kk);
    ua = A\b;
    ua_completo = [0; ua; 0];
    
    % Salva l'errore corrente e i DoFs per il plot finale
    dof_adapt(i) = length(ua); % Numero di DoF interni
    errH1s_adapt_vettore(i) = H1s_1D_error(xnodes_adapt, gradphi, ne, quadrature_type, namegradexactu, ua_completo, alpha);
    
    % 2. ESTIMATE (Calcolo dello stimatore locale eta^2 per ogni elemento)
    eta2 = stimatore_errore_1D(xnodes_adapt, ua_completo, namef, alpha, kk);
    
    % 3. MARK (mark il 25% degli elementi con errore maggiore)
    [~, indici_ordinati] = sort(eta2, 'descend');
    % numero di elementi da marcare (arrotondiamo il 25% fino all'intero più vicino)
    n_marcare = round(0.25 * NS);
    elementi_marcati = indici_ordinati(1:n_marcare);
    
    % 4. REFINE (aggiungiamo i punti medi degli elementi marcati)
    nuovi_nodi = (xnodes_adapt(elementi_marcati) + xnodes_adapt(elementi_marcati+1))/2;
    % uniamo e ordiniamo tutti i nodi eliminando duplicati
    xnodes_adapt = unique([xnodes_adapt; nuovi_nodi]);
end

%% PLOT Grafico Errore in seminorma 1 rispetto ai DoFs

dof_uniforme = zeros(n_mesh, 1);    % DoFs del caso uniforme (punti interni) [xn era la dimensione del sistema]
for idx = 1:n_mesh
    switch kk
        case 1
            dof_uniforme(idx) = nnodes_vector(idx)-2;
        case 2
            dof_uniforme(idx) = (nnodes_vector(idx)-1)*2-1;
    end
end

figure('Name', 'Confronto Uniforme vs Adattativo', 'NumberTitle', 'off');
loglog(dof_uniforme, errH1s_vettore,'Color','#357bb9','Marker','*','LineWidth', 1.5);
hold on;
loglog(dof_adapt, errH1s_adapt_vettore, 'Color','#b9395e','Marker','*', 'LineWidth', 1.5);

% Inseriamo le rette di pendenza teorica come riferimento ottimale
pendenza_ottimale1 = dof_uniforme.^(-1) * (errH1s_adapt_vettore(1) / dof_uniforme(1)^(-1));
loglog(dof_uniforme, pendenza_ottimale1, 'k--', 'LineWidth', 1), hold on;
pendenza_ottimale2 = dof_uniforme.^(-2) * (errH1s_adapt_vettore(1) / dof_uniforme(1)^(-2));
loglog(dof_uniforme, pendenza_ottimale2, 'k:', 'LineWidth', 1);
legend('Mesh Uniforme', 'Mesh Adattativa','Pendenza Teorica N^{-1}', 'Pendenza Teorica N^{-2}', 'Location', 'Best');


title(['Confronto Uniforme vs Adattativo (\alpha = ', num2str(alpha), ', k = ', num2str(kk), ')']);
xlabel('Numero di DoFs');
ylabel('|u-u_h|_{1,\Omega}');
grid on;
axis tight;

%% SALVO DATI PER PLOT
% Per avere i plot questo main va fatto andare 4 volte, facendo tutte le
% possibili combinazioni di alpha e kk

if alpha == 10
    str_alpha = '10';
else
    str_alpha = '43';
end

nome_file = ['dati_confronto_alpha', str_alpha, '_k', num2str(kk), '.mat'];

% Salviamo sia i vettori uniformi che quelli adattativi
save(nome_file, 'dof_uniforme', 'errH1s_vettore', 'dof_adapt', 'errH1s_adapt_vettore','xnodes_adapt', 'alpha', 'kk');
fprintf('Dati salvati correttamente in: %s\n', nome_file);