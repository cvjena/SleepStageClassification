%% Calculates "Fleiss Kappa" between 3 raters (in preprocessed dataset)
load ('E:\Sleep_classification\2_processedData\preprocessed_data.mat', 'daten');    % load preprocessed data

% remove unused data
fieldNamesunsorted = fieldnames(daten);
patterns_to_remove = {'ECG', 'Piezo', 'Sleepstage_3', 'Sleepstage_1', 'Sleepstage_2', 'Sleepstage_4'};
remove_indices = false(length(fieldNamesunsorted), 1);  
for i = 1:length(patterns_to_remove)
    remove_indices = remove_indices | startsWith(fieldNamesunsorted, patterns_to_remove{i});
end
fieldNamesunsorted = fieldNamesunsorted(~remove_indices);  % Only keep entries that do NOT start with the patterns

% sort patients by number
for i = 1:length(fieldNamesunsorted)
    current_field = fieldNamesunsorted{i};  
    last_char = current_field(end);         
    last_digits(i) = str2double(last_char); 
end
[~, sort_index] = sort(last_digits);        
fieldNames = fieldNamesunsorted(sort_index);
Z = [];


for j = 1:3:length(fieldNames)
    if contains(fieldNames{j}, 'Sleepstage_5')      % "Sleepstage_5", "Sleepstage_6" and "Sleepstage_7" contain label information from rater 1,2 & 3
        Rater_1 = daten.(fieldNames{j});            % label rater 1
        Rater_2 = daten.(fieldNames{j+1});          % label rater 2
        Rater_3 =  daten.(fieldNames{j+2});         % label rater 3
    end

    % create matrix (Input for kappa calculation)
    label_of_raters = [Rater_1 Rater_2 Rater_3];
    number_of_label = [];
    for i = 1 : length(label_of_raters)              % over each patient
        zeros = sum(label_of_raters(i, :) == 0);     % # of AS
        ones  = sum(label_of_raters(i, :) == 1);     % # of W
        twos  = sum(label_of_raters(i, :) == -1);    % # of QS
        number_of_label(i,1) = ones;                
        number_of_label(i,2) = zeros;
        number_of_label(i,3) = twos;
    end

    % calculate Fleiss Kappa for echt patient
    [T] = fleiss(number_of_label);
    Z = [Z;T];

end
% SD & SEM of kappa
SD = std(Z(:,1));
meanK = mean(Z(:,1));
SEM = SD./sqrt(length(fieldNames)/3); 

