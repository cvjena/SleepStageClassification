%% Calculate "Cohens Kappa" between manual label and Sleepware label

% load preprocessed data
load('E:\Sleep_classification\2_processedData\preprocessed_data.mat', 'daten');; 

% trim data
pattern = '^(Sleepstage_1|Sleepstage_4)';                           % Regular expression for fields to be kept (only manual and software label)
fields = fieldnames(daten);
keepFields = fields(~cellfun(@isempty, regexp(fields, pattern)));   % Select fields that match the pattern
filteredData = struct();                                            % Create new struct with the selected fields
for i = 1:numel(keepFields)
    filteredData.(keepFields{i}) = data.(keepFields{i});
end

% merge data on all patients
Label_manuell   = [];       
Label_Software  = [];   

for i = 1:2:length(keepFields)
    Label_manuell = [Label_manuell; filteredData.(keepFields{i})];    
    Label_Software = [Label_Software; filteredData.(keepFields{i+1})];
end

% calculate Cohens kappa
kappa = cohensKappa(Label_manuell, Label_Software);
disp(['Kappa between manual and Sleepware label is: ', num2str(kappa)]);
