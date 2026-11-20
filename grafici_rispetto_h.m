%% SCRIPT GRAFICI: SUBPLOT ALPHA (CONFRONTO P1 VS P2 RISPETTO A h)
% Serve aver fatto andare il main con entrambi i valori di alpha e k (4 volte in totale)
clear; clc; close all;

% Controllo esistenza dei file prima di procedere
files_richiesti =  {'dati_confronto_alpha43_k1.mat', 'dati_confronto_alpha43_k2.mat'; ...
                    'dati_confronto_alpha10_k1.mat', 'dati_confronto_alpha10_k2.mat'};
for f = 1:length(files_richiesti)
    if ~exist(files_richiesti{f}, 'file')
        error('Manca il file %s!', files_richiesti{f});
    end
end

%% LOAD dei Dati

% --- CASO ALPHA = 4/3, k = 1 ---
tmp = load('dati_confronto_alpha43_k1.mat');
nnodes_vector   = tmp.dof_uniforme;
err_a43_k1   = tmp.errH1s_vettore;

% --- CASO ALPHA = 4/3, k = 2 ---
tmp = load('dati_confronto_alpha43_k2.mat');
err_a43_k2   = tmp.errH1s_vettore;

% --- CASO ALPHA = 10, k = 1 ---
tmp = load('dati_confronto_alpha10_k1.mat');
err_a10_k1   = tmp.errH1s_vettore;

% --- CASO ALPHA = 10, k = 2 ---
tmp = load('dati_confronto_alpha10_k2.mat');
err_a10_k2   = tmp.errH1s_vettore;

clear tmp; % Puliamo la variabile temporanea dalla memoria


%% DOMINIO e PASSO della mesh
xmin = 0; xmax = 1; 
h = (xmax - xmin) ./ (nnodes_vector - 1);

% Costanti di shift per allineare visivamente le rette teoriche ai dati reali
% (Prendiamo come riferimento l'ordine di grandezza del primo errore)
C_sing_1 = err_a43_k1(1) / h(1);
C_reg_1  = err_a10_k1(1) / h(1);
C_reg_2  = err_a10_k2(1) / (h(1)^2);

%% ==================== COSTRUZIONE GRAFICI ====================
figure('Name', 'Analisi di Convergenza Uniforme', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 480]);

% --- SINISTRA: Alpha = 4/3 (Soluzione Singolare) ---
subplot(1,2,1);
loglog(h, err_a43_k1, 'Color','#357bb9','Marker','*', 'LineWidth', 1.5); hold on;
loglog(h, err_a43_k2, 'Color','#b9395e','Marker','*', 'LineWidth', 1.5);
% Riferimento teorico: O(h) scalato per essere parallelo vicino ai dati
loglog(h, C_sing_1 * h.^1, '--k', 'LineWidth', 1.2);

title('$\alpha = 4/3$ (Soluzione a Regolarit\`a Ridotta)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Passo della mesh $h$', 'Interpreter', 'latex', 'FontSize', 13); 
ylabel('Errore $|u - u_h|_{1,\Omega}$', 'Interpreter', 'latex', 'FontSize', 13);
legend({'Lineari ($P_1$)', 'Quadratici ($P_2$)', 'Rif. $\mathcal{O}(h)$'}, ...
       'Interpreter', 'latex', 'Location', 'Best', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on; axis tight;

% --- DESTRA: Alpha = 10 (Soluzione Regolare) ---
subplot(1,2,2);
loglog(h, err_a10_k1,'Color','#357bb9','Marker','*', 'LineWidth', 1.5); hold on;
loglog(h, err_a10_k2,'Color','#b9395e','Marker','*', 'LineWidth', 1.5);
% Riferimenti teorici: O(h) e O(h^2) scalati
loglog(h, C_reg_1 * h.^1, '--k', 'LineWidth', 1.2);
loglog(h, C_reg_2 * h.^2, ':k', 'LineWidth', 1.2);

title('$\alpha = 10$ (Soluzione Regolare)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Passo della mesh $h$', 'Interpreter', 'latex', 'FontSize', 13); 
ylabel('Errore $|u - u_h|_{1,\Omega}$', 'Interpreter', 'latex', 'FontSize', 13);
legend({'Lineari ($P_1$)', 'Quadratici ($P_2$)', 'Rif. $\mathcal{O}(h)$', 'Rif. $\mathcal{O}(h^2)$'}, ...
       'Interpreter', 'latex', 'Location', 'Best', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on; axis tight;