%% Create arrays from all ECG and piezo data
clear
close all
parent_folder = 'D:\Studie_Ausgewertet_alle\';             % Input: raw data
outputFolder = 'D:\scripts_output\';                       % Output: struct with merged edf files -> "Data_ECG_Piezo.mat"

% add signal processing toolbox and no edf function 
%%
subfolders = dir(parent_folder);
subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));
for folder_idx = 1:length(subfolders)
    current_folder = fullfile(parent_folder, subfolders(folder_idx).name);
    edf_files = dir(fullfile(current_folder, '*].edf')); % choose all edf files
    if isempty(edf_files)
        disp(['Keine .edf-Dateien im Ordner ', current_folder, ' gefunden.']);
        continue;
    end
    disp(['Verarbeite Ordner: ', current_folder]);
    ECG_cellarray = cell(0);
    ECG_cellarrayP = cell(0);
    ECG_array = [];
    Piezo_array = [];
  
    classification_1 = fullfile(current_folder, 'Dataset_A.csv');   % Dataset A
    cv_file_1 = readtable(classification_1);
    cv_file_1 = table2array(cv_file_1);
    classification_2 = fullfile(current_folder, 'Dataset_B.csv');   % Datset B
    cv_file_2 = readtable(classification_2);
    cv_file_2 = table2array(cv_file_2);
    classification_3 = fullfile(current_folder, 'Dataset_C.csv');   % Dataset C
    cv_file_3 = readtable(classification_3);
    cv_file_3 = table2array(cv_file_3);
    fileList = dir(current_folder);
    G3_ind = {};
    for i = 1:length(fileList)
        if contains(fileList(i).name, 'SleepStaging_G3', 'IgnoreCase', true)
            G3_ind{end+1} = fileList(i).name; % Add to matching list
        end
    end
    classification_4 = fullfile(current_folder, G3_ind{1});  % sleepware label
    cv_file_4 = readtable(classification_4);
    cv_file_4 = cv_file_4.Schlafstadium;
    cv_file_4 = replace(cv_file_4, {'WK', 'N2','NS','N', 'REM'}, {'1', '-1', '-1','-1', '0'});
    cv_file_4 = str2double(cv_file_4);
    classification_5 = fullfile(current_folder, 'SleepStaging_CD.csv');  % CD (rater 1)
    cv_file_5 = readtable(classification_5);
    cv_file_5 = cv_file_5.Schlafstadium;
    cv_file_5 = replace(cv_file_5, {'WK', 'N2','NS','N', 'REM'}, {'1', '-1', '-1','-1', '0'});
    cv_file_5 = str2double(cv_file_5);
    classification_6 = fullfile(current_folder, 'SleepStaging_LS.csv');  % LS (rater 2)
    cv_file_6 = readtable(classification_6);
    cv_file_6 = cv_file_6.Schlafstadium;
    cv_file_6 = replace(cv_file_6, {'WK', 'N2','NS','N', 'REM'}, {'1', '-1', '-1','-1', '0'});
    cv_file_6 = str2double(cv_file_6);
    classification_7 = fullfile(current_folder, 'SleepStaging_SJ.csv');  % SJ (rater 3)
    cv_file_7 = readtable(classification_7);
    cv_file_7 = cv_file_7.Schlafstadium;
    cv_file_7 = replace(cv_file_7, {'WK', 'N2','NS','N', 'REM'}, {'1', '-1', '-1','-1', '0'});
    cv_file_7 = str2double(cv_file_7);

    for i = 1:length(edf_files) % length = number of edf-files
        % timetable with only EEG channel
        currentFile = fullfile(current_folder, edf_files(i).name);
        % [~, ECG_signal] = edfread(currentFile, 'targetSignals', 'ECG II'); % Jürgens Version von edfread
        % [~, Piezo_signal] = edfread(currentFile, 'targetSignals', 'Numeric Aux1'); % Jürgens Version von edfread
        ECG_signal = cell2mat(table2array(edfread(currentFile, 'SelectedSignals', 'ECG II')));
        Piezo_signal = cell2mat(table2array(edfread(currentFile, 'SelectedSignals', 'Numeric Aux1')));
        ECG_array = [ECG_array; ECG_signal];
        Piezo_array = [Piezo_array; Piezo_signal];
    end
    % ECG_array = vertcat(ECG_cellarray{:});
    % Piezo_array = vertcat(ECG_cellarrayP{:});
    parts = strsplit(current_folder, '\');
   
    lastPart = parts{end};  % last part of path
    parts = strsplit(lastPart, ' ');
    lastTwoParts = [parts(1),'_', parts(2)];
    lastTwoParts = cell2mat(lastTwoParts);
    part = strsplit(lastTwoParts, '_');
    number = part(2);
    number = cell2mat(number);

    variableNameECG = ['ECG_', number];
    variableNamePiezo = ['Piezo_',  number];
    variableNameSleepstage1 = ['Sleepstage_1_',  number];
    variableNameSleepstage2 = ['Sleepstage_2_',  number];
    variableNameSleepstage3 = ['Sleepstage_3_',  number];
    variableNameSleepstage4 = ['Sleepstage_4_',  number];
    variableNameSleepstage5 = ['Sleepstage_5_',  number];
    variableNameSleepstage6 = ['Sleepstage_6_',  number];
    variableNameSleepstage7 = ['Sleepstage_7_',  number];

    daten.(variableNameECG) = ECG_array;
    daten.(variableNamePiezo) = Piezo_array;
    daten.(variableNameSleepstage1) = cv_file_1;
    daten.(variableNameSleepstage2) = cv_file_2;
    daten.(variableNameSleepstage3) = cv_file_3;
    daten.(variableNameSleepstage4) = cv_file_4;      % sleepware label
    daten.(variableNameSleepstage5) = cv_file_5;
    daten.(variableNameSleepstage6) = cv_file_6;
    daten.(variableNameSleepstage7) = cv_file_7;
    
end
filename3 = [outputFolder, 'Data_ECG_Piezo.mat']; % struct with Patient1:ECG,Piezo, Sleepstage; Patient2...
save(filename3, 'daten');
disp(['daten were saved in ', filename3]);
%% number of epochs

totalSize = 0;
fields = fieldnames(daten);
for i = 1:length(fields)
    fieldName = fields{i};
    if contains(fieldName, 'Sleepstage_1')
        totalSize = totalSize + numel(daten.(fieldName));
    end
end

disp(['number of epochs is : ', num2str(totalSize)]);
