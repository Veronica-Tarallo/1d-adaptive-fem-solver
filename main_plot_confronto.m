%% GRAFICI ADATTATIVO VS UNIFORME
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

% Lista dei file da caricare
%   Riga 1: Alpha = 4/3, k=1 e k=2
%   Riga 2: Alpha = 10, k=1 e k=2
files = {...
    'dati_confronto_alpha43_k1.mat', 'dati_confronto_alpha43_k2.mat'; ...
    'dati_confronto_alpha10_k1.mat', 'dati_confronto_alpha10_k2.mat'};
 
titoli = {...
    '$\alpha = 4/3$, $P_1$ (Lineare)', '$\alpha = 4/3$, $P_2$ (Quadratico)'; ...
    '$\alpha = 10$, $P_1$ (Lineare)', '$\alpha = 10$, $P_2$ (Quadratico)'};
fig_all = figure('Name', 'Confronto Globale Adattativo vs Uniforme');

for r = 1:2
    for c = 1:2
        % Carica i dati correnti
        data = load(files{r,c});
        
        % Seleziona il quadrante corretto del subplot
        idx_sub = (r-1)*2 + c;
        ax = subplot(2, 2, idx_sub);
        
        % Grafico dati Uniformi
        loglog(data.dof_uniforme, data.errH1s_vettore, 'Color', '#357bb9', 'Marker', '*', 'LineWidth', 1.5);
        hold on;
        % Grafico dati Adattativi
        loglog(data.dof_adapt, data.errH1s_adapt_vettore, 'Color', '#b9395e', 'Marker', '*', 'LineWidth', 1.5);
        
        % Calcolo e plot delle pendenze teoriche locali
        pendenza1 = data.dof_uniforme.^(-1) * (data.errH1s_vettore(1) / data.dof_uniforme(1)^(-1));
        pendenza2 = data.dof_uniforme.^(-2) * (data.errH1s_vettore(1) / data.dof_uniforme(1)^(-2));
        
        loglog(data.dof_uniforme, pendenza1, '--k', 'LineWidth', 1.2);
        loglog(data.dof_uniforme, pendenza2, ':k', 'LineWidth', 1.2);
        
        set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
        title(titoli{r,c}, 'Interpreter', 'latex', 'FontSize', 14);
        xlabel('Gradi di Libert\`a (DoFs)', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$|u - u_h|_{1,\Omega}$', 'Interpreter', 'latex', 'FontSize', 12);
        
        grid on;
        axis tight;
        
        % Inseriamo la legenda solo nel primo quadrante
        if idx_sub == 1
            lgd = legend('Mesh Uniforme', 'Mesh Adattativa', ...
                         'Pendenza $\mathcal{O}(N^{-1})$', 'Pendenza $\mathcal{O}(N^{-2})$', ...
                         'Location', 'best');
            set(lgd, 'Interpreter', 'latex', 'FontSize', 10);
        end
    end
end