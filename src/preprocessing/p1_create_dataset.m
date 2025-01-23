%% Create dataset A, B and C (consensus label variants) 
% Plot: 3 different label of sleep stages


% import add-ons: Fleiss
parent_folder = 'C:\Users\user\Downloads\Studie_Ausgewertet_volstaendig\Studie_Ausgewertet_alle\';       % Input folder containing folders for each patient with .csv files for each rater

%%
subfolders = dir(parent_folder);
subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));

Z = [];
mittel = 0;
mittelnum_W = 0;
mittelnum_AS = 0;
mittelnum_QS = 0;
sd = 0;
sj = 0;
dj = 0;
fs = 0;
fd = 0;
fj = 0;
kappa_patient = [];
num_epochs = 0;
num_NaN = 0;
num12 = 0;
num13 = 0;
num23 = 0;
label_r1 = [];
label_r2 = [];
label_r3 = [];
num_equal_label = 0;

for folder_idx = 1:length(subfolders)
    % read data
    current_folder = fullfile(parent_folder, subfolders(folder_idx).name);
    titelname = subfolders(folder_idx).name(1:9);
    titelname = strrep(titelname, 'Studie', 'Patient');
    csv_11 = dir(fullfile(current_folder, '*_LS.csv'));
    csv_1 = [current_folder, '\',csv_11.name];
    csv_1 = readtable(csv_1);
    csv_1 = csv_1.Schlafstadium;
    csv_1 = replace(csv_1, {'WK', 'N2','N', 'REM'}, {'1', '-1', '-1', '0'});
    csv_1 = str2double(csv_1);
    num_epochs=num_epochs+length(csv_1);
    
    csv_22 = dir(fullfile(current_folder, '*_CD.csv'));
    csv_2 = [current_folder, '\',csv_22.name];
    csv_2 = readtable(csv_2);
    csv_2 = csv_2.Schlafstadium;
    csv_2 = replace(csv_2, {'WK', 'N2','N', 'REM'}, {'1', '-1', '-1', '0'});
    csv_2 = str2double(csv_2);
    
    csv_33 = dir(fullfile(current_folder, '*_SJ.csv'));    
    csv_3 = [current_folder, '\',csv_33.name];    
    csv_3 = readtable(csv_3);    
    csv_3 = csv_3.Schlafstadium;    
    csv_3 = replace(csv_3, {'WK', 'N2','N', 'REM'}, {'1', '-1', '-1', '0'});
    csv_3 = str2double(csv_3);
   
    csv_4 = NaN(length(csv_1),1);       
    csv_5 = NaN(length(csv_1),1);
    csv_6 = NaN(length(csv_1),1);
    
    % create dataset (consensus label)
    for i=1:length(csv_1)
        if csv_1(i) == csv_2(i) && csv_1(i) == csv_3(i)
            csv_4(i) = csv_1(i);
            csv_5(i) = csv_1(i);
            csv_6(i) = csv_1(i);
            num_equal_label = num_equal_label+1;
        elseif csv_1(i)==csv_2(i)
            csv_4(i) = csv_1(i);
            num12 = num12+1;
        elseif csv_1(i)==csv_3(i)
            csv_4(i) = csv_1(i);
            num13 = num13+1;
        elseif csv_2(i)==csv_3(i)
            csv_4(i)=csv_2(i);
            num23 = num23+1;          
        else
            csv_4(i) = NaN;
            csv_5(i) = NaN;
            num_NaN = num_NaN+1;            
        end
    end
    for i=1:length(csv_1)
        if csv_1(i) == 0 | csv_2(i)==0 | csv_3(i) ==0
           csv_6(i)=0;            
        end
    end

    % save
    pfad = [current_folder, '\Dataset_A.csv'];  % Dataset A
    writematrix(csv_5, pfad);
    pfad = [current_folder, '\Dataset_B.csv'];  % Dataset B
    writematrix(csv_4, pfad);   
    pfad = [current_folder, '\Dataset_C.csv'];  % Dataset C
    writematrix(csv_6, pfad);


    % compare rater (check whether a rater stands out particularly
    % strongly)
    % rater = D., J., S.
    sd = sd + sum(csv_1==csv_2);
    sj = sj + sum(csv_1==csv_3);
    dj = dj + sum(csv_2==csv_3);

    fs = fs + sum(csv_2==csv_3 & csv_1~=csv_2);
    fd = fd + sum(csv_1==csv_3& csv_1~=csv_2);
    fj = fj + sum(csv_2==csv_1& csv_1~=csv_3);

    % create matrix for kappa calculation
    x = [csv_1 csv_2 csv_3];
    matrix = [];
    for i = 1 : length(x)
        num_zeros = sum(x(i, :) == 0);
        num_ones = sum(x(i, :) == 1);
        num_twos = sum(x(i, :) == -1);
        matrix(i,1) = num_ones;
        matrix(i,2) = num_zeros;
        matrix(i,3) = num_twos;
    end
    
    % calculate kappa (for all epochs -  not trimmed data)
    T = fleiss(matrix);
    mittel = mittel + T.Fleiss_k;
    Z = [Z;T];
  
    kappa_patient(folder_idx).name = titelname;
    kappa_patient(folder_idx).kappa = T.Fleiss_k;
  
    label_r1 = [label_r1; csv_1];
    label_r2 = [label_r2; csv_2];
    label_r3 = [label_r3; csv_3];
end


meanKappa = mittel/length(subfolders);
meanKappanum_W = mittelnum_W/length(subfolders);
meanKappanum_AS = mittelnum_AS/length(subfolders);
meanKappanum_QS = mittelnum_QS/length(subfolders);
percentage_NaN = num_NaN/num_epochs;                       % numer of epochs without label
num2 = num12 + num23 + num13;
output = ['mean Kappa: ',num2str(meanKappa)];
disp(output);


% SD & SEM of kappa
SD = std(Z(:,1));
meanK = mean(Z(:,1));
SEM = SD./sqrt(28);

%% Figure: compare rater 1,2,3

figure('WindowState', 'maximized');
tiledlayout(6,1)

label_r1 = csv_3;
label_r2 = csv_2;
label_r3 = csv_1;

% define colours
colors = [1, 0, -1];
bar_colors = [215/255, 104/255, 45/255; 
              89/255, 189/255, 107/255;
              119/255, 206/255, 219/255];

% legend
labels = {'W', 'AS', 'QS'};

% figure
ax1 = nexttile;
hold on;

% bars for legend
legend_entries = [];
legend_labels = {};

% draw
draw_rectangle = @(x, y, width, height, color) rectangle('Position', [x, y, width, height], 'FaceColor', color, 'EdgeColor', 'none');


for i = 1:length(label_r1)
    value = label_r1(i);
    if isnan(value)
        continue;
    end
    color_idx = find(colors == value);
    draw_rectangle(i-0.5, 0, 1, 0.4, bar_colors(color_idx, :));
    
    if isempty(legend_labels) || ~any(strcmp(legend_labels, labels{color_idx}))
        legend_entries = [legend_entries, bar(i, NaN, 'FaceColor', bar_colors(color_idx, :), 'EdgeColor', 'none')];
        legend_labels = [legend_labels, labels{color_idx}];
    end
end

for i = 1:length(label_r2)
    value = label_r2(i);
    if isnan(value)
        continue;
    end
    color_idx = find(colors == value);
    draw_rectangle(i-0.5, 0.45, 1, 0.4, bar_colors(color_idx, :));
    
    if isempty(legend_labels) || ~any(strcmp(legend_labels, labels{color_idx}))
        legend_entries = [legend_entries, bar(i, NaN, 'FaceColor', bar_colors(color_idx, :), 'EdgeColor', 'none')];
        legend_labels = [legend_labels, labels{color_idx}];
    end
end

for i = 1:length(label_r3)
    value = label_r3(i);
    if isnan(value)
        continue; 
    end
    color_idx = find(colors == value);
    draw_rectangle(i-0.5, 0.9, 1, 0.4, bar_colors(color_idx, :));
    
    if isempty(legend_labels) || ~any(strcmp(legend_labels, labels{color_idx}))
        legend_entries = [legend_entries, bar(i, NaN, 'FaceColor', bar_colors(color_idx, :), 'EdgeColor', 'none')];
        legend_labels = [legend_labels, labels{color_idx}];
    end
end

title('Compare rater annotation');
ylim([0, 1.3]); 
xlim([0, length(label_r1) + 1]); 
set(gca, 'xtick', []); 
set(gca, 'ytick', []);

% y label
text(-1, 0.23, 'Rater 3', 'Rotation', 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontSize', 25);
text(-1, 0.66, 'Rater 2', 'Rotation', 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontSize', 25);
text(-1, 1.08, 'Rater 1', 'Rotation', 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontSize', 25);

% add legend
legend(legend_entries, legend_labels, 'Location', 'eastoutside','FontSize',25);
legend boxoff;
hold off;

%% Figure: Compare dataset A, B, C

figure('WindowState', 'maximized');
tiledlayout(6,1)

label_r1 = csv_6;
label_r2 = csv_4;
label_r3 = csv_5;

% define colour
colors = [1, 0, -1];
bar_colors = [215/255, 104/255, 45/255; 
              89/255, 189/255, 107/255;
              119/255, 206/255, 219/255]; 

% legend
labels = {'W', 'AS', 'QS'};

ax1 = nexttile;
hold on;

legend_entries = [];
legend_labels = {};

draw_rectangle = @(x, y, width, height, color) rectangle('Position', [x, y, width, height], 'FaceColor', color, 'EdgeColor', 'none');

for i = 1:length(label_r1)
    value = label_r1(i);
    if isnan(value)
        continue; 
    end
    color_idx = find(colors == value);
    draw_rectangle(i-0.5, 0, 1, 0.4, bar_colors(color_idx, :));
    
    if isempty(legend_labels) || ~any(strcmp(legend_labels, labels{color_idx}))
        legend_entries = [legend_entries, bar(i, NaN, 'FaceColor', bar_colors(color_idx, :), 'EdgeColor', 'none')];
        legend_labels = [legend_labels, labels{color_idx}];
    end
end


for i = 1:length(label_r2)
    value = label_r2(i);
    if isnan(value)
        continue; 
    end
    color_idx = find(colors == value);
    draw_rectangle(i-0.5, 0.45, 1, 0.4, bar_colors(color_idx, :));
    
    if isempty(legend_labels) || ~any(strcmp(legend_labels, labels{color_idx}))
        legend_entries = [legend_entries, bar(i, NaN, 'FaceColor', bar_colors(color_idx, :), 'EdgeColor', 'none')];
        legend_labels = [legend_labels, labels{color_idx}];
    end
end

for i = 1:length(label_r3)
    value = label_r3(i);
    if isnan(value)
        continue; 
    end
    color_idx = find(colors == value);
    draw_rectangle(i-0.5, 0.9, 1, 0.4, bar_colors(color_idx, :));
    
    if isempty(legend_labels) || ~any(strcmp(legend_labels, labels{color_idx}))
        legend_entries = [legend_entries, bar(i, NaN, 'FaceColor', bar_colors(color_idx, :), 'EdgeColor', 'none')];
        legend_labels = [legend_labels, labels{color_idx}];
    end
end

title('Compare dataset');
ylim([0, 1.3]);
xlim([0, length(label_r1) + 1]); 
set(gca, 'xtick', []); 
set(gca, 'ytick', []);

text(-1, 0.24, 'Dataset C', 'Rotation', 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontSize', 25);
text(-1, 0.64, 'Dataset B', 'Rotation', 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontSize', 25);
text(-1, 1.08, 'Dataset A', 'Rotation', 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontSize', 25);

[~, sort_idx] = ismember(labels, legend_labels);
[~, sorted_idx] = sort(sort_idx);
sorted_entries = legend_entries(sorted_idx);
sorted_labels = legend_labels(sorted_idx);

% add legend
legend(sorted_entries, sorted_labels, 'Location', 'eastoutside', 'FontSize',25);
legend boxoff;
hold off;



%% % sleep stage (all data)

subfolders = dir(parent_folder);
subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));

num_epochs = 0;
label_r1 = [];
label_r2 = [];
label_r3 = [];
percentage_W = [];
percentage_AS = [];
percentage_QS = [];

for folder_idx = 1:length(subfolders)
    current_folder = fullfile(parent_folder, subfolders(folder_idx).name);

    csv_11 = dir(fullfile(current_folder, '*_LS.csv'));
    csv_1 = [current_folder, '\',csv_11.name];
    csv_1 = readtable(csv_1);
    csv_1 = csv_1.Schlafstadium;
    csv_1 = replace(csv_1, {'WK', 'N2','N', 'REM'}, {'1', '-1', '-1', '0'});
    csv_1 = str2double(csv_1);
    num_epochs = num_epochs + length(csv_1);
    
    csv_22 = dir(fullfile(current_folder, '*_CD.csv'));
    csv_2 = [current_folder, '\',csv_22.name];
    csv_2 = readtable(csv_2);
    csv_2 = csv_2.Schlafstadium;
    csv_2 = replace(csv_2, {'WK', 'N2','N', 'REM'}, {'1', '-1', '-1', '0'});
    csv_2 = str2double(csv_2);
    
    csv_33 = dir(fullfile(current_folder, '*_SJ.csv'));    
    csv_3 = [current_folder, '\',csv_33.name];    
    csv_3 = readtable(csv_3);    
    csv_3 = csv_3.Schlafstadium;    
    csv_3 = replace(csv_3, {'WK', 'N2','N', 'REM'}, {'1', '-1', '-1', '0'});
    csv_3 = str2double(csv_3);  
    label_r1 = [label_r1; csv_1];
    label_r2 = [label_r2; csv_2];
    label_r3 = [label_r3; csv_3];
    label_r123 = 0;
    label_r123 = [label_r1; label_r2; label_r3];
    num_W = sum(label_r123 == 1);
    num_AS = sum(label_r123 == 0);
    num_QS = sum(label_r123 == -1);

    percentage_W = [percentage_W,num_W];
    percentage_AS = [percentage_AS,num_AS];
    percentage_QS = [percentage_QS,num_QS];
    
end
label_r123 = [label_r1; label_r2; label_r3];
num_W = sum(label_r123 == 1);
num_AS = sum(label_r123 == 0);
num_QS = sum(label_r123 == -1);

percentage_W = (num_W*100)/(num_epochs*3);      % (3 rater)
percentage_AS = (num_AS*100)/(num_epochs*3);
percentage_QS = (num_QS*100)/(num_epochs*3);

%% Figure: pie chart % 3,2,1,0 rater agreement (all data)

daten = [67.914, 9.046, 13.511, 8.326, 1.296];

labels = {'3 Agreements', '2 Agreements (1 & 2)', '2 Agreements (1 & 3)', '2 Agreements (2 & 3)', '0 Agreements'};

% define colours
colours_stage = [0, 0.447, 0.741;
            0.226, 0.572, 0.800;
            0.452, 0.697, 0.858;
            0.678, 0.822, 0.917;
            0.900, 0.945, 0.974];

figure; 
h = pie(daten);

patchHandles = findobj(h, 'Type', 'Patch');
for i = 1:length(patchHandles)
    set(patchHandles(i), 'FaceColor', colours_stage(i, :));
end

lgd = legend(labels, 'Location', 'eastoutside');
legend boxoff;

set(lgd, 'FontSize', 28); 

textHandles = findobj(h, 'Type', 'Text');
set(textHandles, 'FontSize', 28); 
