clc
clear all
close all 
global dates
dates = datetime(2020,3,1):days(1):datetime(2023,6,30);


plot_scen_data('output/sim_diff/vac_1217_cnt0_0_0_0_100_n75_s3_belgium/cumul_deaths_incr0_0_0_0_100_vac_vaccine_uptake_extrabooster_vscenbelgium_fall23campaign_0flu_biv_belgium.mat', ...
               'plot_data/deaths.csv',true);
xlabel('Date')
ylabel('Cumulative deaths')
fig=gcf();
fig.UserData.filename = 'vac_deaths';
fig.UserData.orientation = 'landscape';

plot_scen_data('output/sim_diff/vac_1217_cnt0_0_0_0_100_n75_s3_belgium/cumul_cases_incr0_0_0_0_100_vac_vaccine_uptake_extrabooster_vscenbelgium_fall23campaign_0flu_biv_belgium.mat', ...
               'plot_data/cases.csv',true);
xlabel('Date')
ylabel('Cumulative cases')
fig=gcf();
fig.UserData.filename = 'vac_cases';
fig.UserData.orientation = 'landscape';

figure;
plot_delta_scen('delta_cumul_deaths',65,1)
xlabel('Date')
ylabel('Averted deaths')
plot_delta_scen('delta_cumul_deaths',800,2)
xlabel('Date')
ylabel('Averted deaths')
fig=gcf();
fig.UserData.filename = 'recap_deaths';

plot_scen_data('output/sim_diff/vac_1217_cnt0_0_0_0_100_n75_s3_belgium/hosp_adm_cum0_0_0_0_100_vac_vaccine_uptake_extrabooster_vscenbelgium_fall23campaign_0flu_biv_belgium', ...
   'plot_data/new_hosp2.csv',true)
xlabel('Date')
ylabel('Cumulative hospital admissions')
fig=gcf();
fig.UserData.filename = 'vac_new_hosp';
fig.UserData.orientation = 'landscape';
plot_scen_data('output/sim_diff/vac_1217_cnt0_0_0_0_100_n75_s3_belgium/hosp_load_incr0_0_0_0_100_vac_vaccine_uptake_extrabooster_vscenbelgium_fall23campaign_0flu_biv_belgium', ...
   'plot_data/hosp_load.csv',false)
xlabel('Date')
ylabel('Hospital load')
fig=gcf();
fig.UserData.filename = 'vac_hospload';
fig.UserData.orientation = 'landscape';
plot_scen_data('output/sim_diff/vac_1217_cnt0_0_0_0_100_n75_s3_belgium/icu_load_incr0_0_0_0_100_vac_vaccine_uptake_extrabooster_vscenbelgium_fall23campaign_0flu_biv_belgium', ...
   'plot_data/ICU_load.csv',false)
xlabel('Date')
ylabel('ICU load')
fig=gcf();
fig.UserData.filename = 'vac_icuload';
fig.UserData.orientation = 'landscape';

figure;
plot_delta_scen('delta_hosp_new_cumul',65,1)
xlabel('Date')
ylabel('Absolute difference')
plot_delta_scen('delta_hosp_new_cumul',800,2)
xlabel('Date')
ylabel('Absolute difference')
fig=gcf();
fig.UserData.filename = 'recap_newhosp';

figure;
plot_delta_scen('delta_hosp_load',65,1)
xlabel('Date')
ylabel('Absolute difference')
plot_delta_scen('delta_hosp_load',800,2)
xlabel('Date')
ylabel('Absolute difference')
fig=gcf();
fig.UserData.filename = 'recap_hospload';

figure;
plot_delta_scen('delta_icu_load',65,1)
xlabel('Date')
ylabel('Absolute difference')
plot_delta_scen('delta_icu_load',800,2)
xlabel('Date')
ylabel('Absolute difference')
fig=gcf();
fig.UserData.filename = 'recap_icu';

figure;
plot_delta_scen('delta_icu_load',65,1)
xlabel('Date')
ylabel('Absolute difference')
ylim([-1000 3000])
plot_delta_scen('delta_icu_load',800,2)
xlabel('Date')
ylabel('Absolute difference')
ylim([-1000 3000])
fig=gcf();
fig.UserData.filename = 'recap_icu2';

figure;
plot_delta_scen('delta_cumul_cases',65,1)
xlabel('Date')
ylabel('Averted cases')
ylim([-8e6 8e6]);
plot_delta_scen('delta_cumul_cases',800,2)
xlabel('Date')
ylabel('Averted cases')
ylim([-8e6 8e6]);
fig=gcf();
fig.UserData.filename = 'recap_case';

figure;
plot_delta_scen('social_coef',65,1)
xlabel('Date')
ylabel('s_{bar} (t)')
ylim([0 1])
plot_delta_scen('social_coef',800,2)
xlabel('Date')
ylabel('s_{bar} (t)')
ylim([0 1])
fig=gcf();
fig.UserData.filename = 'recap_social';
% Définition des paramètres


base_path = 'output/sim_diff/';
tag = 'social_coef';
stitles={'Original barometer','Less stringent barometer'};
options = [65, 800]; % Deux valeurs de 'option'
ab = [0 0; 0.1 0.3; 0.2 0.6; 0.4 0.6; 0.2 0.9; 0.4 0.9]; % 6 paires de (a, b)
titles = { 'no barometer', ...
           'orange=0.3, red=0.1', ...
           'orange=0.6, red=0.2', ...
           'orange=0.6, red=0.4', ...
           'orange=0.9, red=0.2', ...
           'orange=0.9, red=0.4' };

% Création de la figure principale
figure;
mainLayout = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
  % 3 lignes, 2 colonnes pour chaque option
%sut(2) = tiledlayout( 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact','Parent',mainLayout);  % 3 lignes, 2 colonnes pour chaque option


% Itération sur les options
for opt_idx = 1:length(options)
    option = options(opt_idx);

    % Recherche des fichiers correspondant à l'option
    search_pattern = sprintf('**/%s*%d*.mat', tag, option);
    files = dir(fullfile(base_path, search_pattern));

    % Recherche des fichiers incluant "nobarometer"
    search_pattern_nobar = sprintf('**/%s*nobarometer*.mat', tag);
    files_nobar = dir(fullfile(base_path, search_pattern_nobar));

    % Combinaison des fichiers
    files = [files_nobar; files];
    
    % Construction des chemins des fichiers
    file_paths = fullfile({files.folder}, {files.name});
    
    % Affichage des fichiers trouvés
    disp(['Found files for option ', num2str(option), ':']);
    disp(file_paths);

    if isempty(file_paths)
        warning(['No files found for option ', num2str(option)]);
        continue;
    end

    % Chargement des données
    data_structs = cellfun(@load, file_paths, 'UniformOutput', false);
    data_matrices = cellfun(@(s) struct2cell(s), data_structs, 'UniformOutput', false);
    data_matrices = cellfun(@(c) c{1}, data_matrices, 'UniformOutput', false); % Extraction de la première variable

    % Création du sous-layout (3x2) pour chaque option (colonne)
    sut=tiledlayout( 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact','Parent',mainLayout);
    sut.Layout.Tile=opt_idx;
    sut.Title.String=stitles{opt_idx};
    sut.Title.FontSize=24;
    sut.Title.FontWeight='bold';
    % Génération des graphiques (1 par couple (a, b))
    for i = 1:min(6, length(data_matrices))
        coef = data_matrices{i};
        total_elements = size(coef, 2);

        % Définition des bornes dynamiques
        a = ab(i,1);
        b = ab(i,2);

        % Calcul des proportions
        percent_below_a = sum(coef <= a, 2) / total_elements * 100;
        percent_between_a_b = sum(coef > a & coef <= b, 2) / total_elements * 100;
        percent_above_b = sum(coef > b, 2) / total_elements * 100;

        % Construction des données
        data = [percent_below_a, percent_between_a_b, percent_above_b];

        % Création du sous-graphique
        ax = nexttile(sut);
        h = bar(dates, data, 'stacked','EdgeAlpha',0,'BarWidth', 1);
        
        % Application des couleurs
        h(1).FaceColor = [1 0 0];      % Rouge pour ≤ a
        h(2).FaceColor = [1 0.5 0];    % Orange pour a < x ≤ b
        h(3).FaceColor = [1 0.95 0.3]; % Jaune plus terne pour > b
        h(1).EdgeColor='none';
        h(2).EdgeColor='none';
        h(3).EdgeColor='none';
        xlim([datetime(2021,1,1), datetime(2023,7,1)]);
        ylim([0 100]);
        grid off
        
        % Ajustement des axes pour éviter la redondance
        if (mod(i, 2) ~= 1 || opt_idx==2) % Cacher l'axe Y des colonnes 2 et 3
            ax.YAxis.Visible = 'off';
        else
            ylabel('%')
        end
        if i <= 4 % Cacher l'axe X des lignes du haut
            ax.XAxis.Visible = 'off';
        else
            xlabel('Date')
        end

        % Titre du sous-graphique
        title({titles{i};''}, 'FontSize', 16);
    end
    
end
fig = gcf();
fig.UserData.filename = 'probacode';



% Définition du dossier d'exportation
export_folder = 'exported_figures';
if ~exist(export_folder, 'dir')
    mkdir(export_folder);
end

% Parcours des figures ouvertes
figHandles = findall(0, 'Type', 'figure');

for i = 1:length(figHandles)
    fig = figHandles(i);
    ax = findall(fig, 'Type', 'axes');
    for j = 1:length(ax)
    ax(j).YAxis.TickLabelFormat = '%g';
    end
    % Récupérer les infos stockées dans UserData
    if isfield(fig.UserData, 'filename')
        filename = fig.UserData.filename;
    else
        filename = sprintf('figure_%02d', i); % Nom par défaut
    end
    
    if isfield(fig.UserData, 'orientation')
        orientation = fig.UserData.orientation;
    else
        orientation = 'landscape'; % Orientation par défaut
    end

    % Définition du format papier
    if strcmp(orientation, 'portrait')
        paperSize = [21, 29.7]; % A4 Portrait (cm)
    else
        paperSize = [18, 12]; % A4 Paysage (cm)
        if(i==1) 
            paperSize=[60,40];
        end
    end
    
    % Ajustement de la taille de la figure
    margin = 0; % Marge en cm
    figSize = paperSize - margin;
    fig.Units = 'centimeters';
    fig.Position = [1, 1, figSize(1), figSize(2)];
    
    % Paramètres d'impression
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = paperSize;
    fig.PaperPosition = [-0.2, 0, paperSize(1)+1.25, paperSize(2)];

    % Amélioration de la lisibilité
    ax = findall(fig, 'Type', 'axes');
    if (i~=1)
    for j = 1:length(ax)
        ax(j).FontSize = 7;
        ax(j).LabelFontSizeMultiplier = 1;
        %ax(j).Position(4) = ax(j).Position(4) * 1.1;
        ax(j).XAxis.TickLabelFormat = 'MM-yy';
    end
    else
    for j = 1:length(ax)
        ax(j).FontSize = 20;
        ax(j).LabelFontSizeMultiplier = 1.5;
        ax(j).XAxis.TickLabelFormat = 'MM-yy';
    end
    end

    
    % Exportation
    filepath = fullfile(export_folder, sprintf('%s.pdf', filename));
    print(fig, filepath, '-dpdf');
    
    fprintf('✅ Figure exportée : %s (%s)\n', filename, orientation);
end

disp('✅ Exportation terminée avec succès !');


function plot_scen_data(scen_path, data_path,cum_tag)
    global dates
    % Load scenario data
    scen_mort = load(scen_path);
    scen_mort_agg = [median(scen_mort.data, 2), ...
                     quantile(scen_mort.data, 0.025, 2), ...
                     quantile(scen_mort.data, 0.975, 2)];
    
    % Generate date axis matching the scenario data
    % Load real-world data
    data = readtable(data_path, 'TextType', 'string', 'Delimiter', 'comma');

    % Liste des mois en français et en anglais
months_fr = ["janv.", "févr.", "mars", "avr.", "mai", "juin", "juil.", ...
             "août", "sept.", "oct.", "nov.", "déc."];
months_en = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", ...
             "Aug", "Sep", "Oct", "Nov", "Dec"];

% Remplacer les noms de mois en français par leur équivalent en anglais
for i = 1:length(months_fr)
    data{:,1} = replace(data{:,1}, months_fr(i), months_en(i));
end

% Convertir les dates en format datetime avec le bon format
data{:,1} = datetime(data{:,1}, 'InputFormat', 'd MMM yyyy', 'Locale', 'en_US');

    % Plot the median
    figure; hold on;
    plot(dates, scen_mort_agg(:,1), 'r', 'LineWidth', 1.5);

    % Plot the confidence interval
    fill([dates, fliplr(dates)], [scen_mort_agg(:,3)', fliplr(scen_mort_agg(:,2)')], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot the real data
    if (cum_tag==true)
    plot(datetime(data{:,1}), cumsum(data{:, 2}), 'ok','MarkerSize',8);
    else
    plot(datetime(data{:,1}), data{:, 2}, 'ok','MarkerSize',8);
    end
    % Set x-axis limits
    xlim([datetime(2021,1,1), datetime(2023,7,1)]);
end


function plot_delta_scen(tag, option,subopt)
    global dates
    legend_list = { 'no barometer', ...
                'orange=0.3, red=0.1', ...
                'orange=0.6, red=0.2', ...
                'orange=0.6, red=0.4', ...
                'orange=0.9, red=0.2', ...
                'orange=0.9, red=0.4' };

    % Validate input
    if nargin < 2
        error('Usage: plot_delta_scen(tag, option), where option is 65 or 800.');
    end
    
    % Search recursively for matching files
    base_path = 'output/sim_diff/';
    search_pattern = sprintf('**/%s*%d*.mat', tag, option);
    files = dir(fullfile(base_path, search_pattern));
    
    % Search for files including "nobarometer"
    search_pattern_nobar = sprintf('**/%s*nobarometer*.mat', tag);
    files_nobar = dir(fullfile(base_path, search_pattern_nobar));
    
    % Combine file lists
    files = [files_nobar; files];
    
    % Build the full list of file paths
    file_paths = fullfile({files.folder}, {files.name});
    
    % Display found files
    disp('Found files:');
    disp(file_paths);
    
    if isempty(file_paths)
        warning('No files found for the given tag and option.');
        return;
    end
    
    % Load the data from the files
    data_structs = cellfun(@load, file_paths, 'UniformOutput', false);
    
    % Extract matrices assuming same field structure
    data_matrices = cellfun(@(s) struct2cell(s), data_structs, 'UniformOutput', false);
    data_matrices = cellfun(@(c) c{1}, data_matrices, 'UniformOutput', false); % Extract first variable
    
    % Compute statistics
    medians = cellfun(@(M) median(M, 2), data_matrices, 'UniformOutput', false);
    quantiles = cellfun(@(M) quantile(M, [0.025, 0.975], 2), data_matrices, 'UniformOutput', false);
    
    % Plot
    subplot(1,2,subopt)
    hold on;
    colors = lines(7); % Define colors
    num_series = length(medians);
    ax = gca;
    if (subopt==1)
    title(ax, 'a) Original barometer', 'Units', 'normalized', 'Position', [0, 1.05, 0],'HorizontalAlignment','left');
    else
    title(ax, 'b) Less stringent barometer', 'Units', 'normalized', 'Position', [0, 1.05, 0],'HorizontalAlignment','left');
    end
    
    for i = 1:num_series
        dates_col = dates(:);
        lower_bound = quantiles{i}(:,1);
        upper_bound = quantiles{i}(:,2);
        color = colors(i+1,:);
        
        % Plot median curve
        plot(dates_col, medians{i}, 'Color', color, 'LineWidth', 1.5);
        
        % Plot uncertainty area without adding it to legend
        fill([dates_col; flipud(dates_col)], [lower_bound; flipud(upper_bound)], ...
             color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    
    xlim([datetime(2021,1,1), datetime(2023,7,1)]);
if(subopt==1)
legend(legend_list,'Location','southoutside');
lgd=legend;
lgd.Position = [0.1 0.02 0.8 0.03]; % Légende sous le graphique
lgd.Orientation = 'horizontal'; % Légende alignée horizontalement
lgd.NumColumns = 6; % Répartir la légende sur plusieurs colonnes
lgd.FontSize = 5; % Réduire la taille de la police
lgd.Box = 'off'; % Supprimer le cadre de la légende pour plus de compacité
end
    hold off;
end