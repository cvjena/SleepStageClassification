%% Calculate Accuracy, Sensitivity, Specificity, and Kappa per Patient

% load SVM model --> Workspace
addpath('E:\SleepClassification\scripts\Workspace\SVM\Dataset A\');
load('ECG+Piezo_Modell.mat');

% load feature
inputFolder = 'D:\Sleep_classification\3_feature\Combination';
fileList = dir(fullfile(inputFolder, '*.csv'));

% Calculate proporton of predicted and observed sleep states
name = [];
Acc_Pat = [];
Sens_Pat = [];
Spez_Pat = [];
Kappa_Pat = [];

for fileIdx = 1:length(fileList)
    % read feature table
    currentFile = fullfile(inputFolder, fileList(fileIdx).name);
    addpath('C:\Users\demme\MATLAB Drive\SleepClassification\scripts\Sleep_analysis\');
    [predicted_label, manual_label] = predict_epochs(currentFile, SVMModel);   % Fkt wendet SVM Model auf Daten an


    kappa_p = cohensKappa(predicted_label, manual_label);    % use "simple Cohens kappa"  

    indexToKeepW = manual_label == 1;
    indexToKeepNW = manual_label ~= 1;
    TP = sum(predicted_label(indexToKeepW) == 1);                               % wach in T und P
    TN = sum(predicted_label(indexToKeepNW) ~= 1);                              % nicht W in T und nicht W in P
    FP = sum(predicted_label(indexToKeepNW) == 1);                              % nicht wach in T aber W in P
    FN = sum(predicted_label(indexToKeepW) ~= 1);
    ACCW = (TP+TN)/(TP+TN+FP+FN);
    SensW = TP/(TP+FN);
    SpesW = TN/(TN+FP);

    indexToKeepAS = manual_label == 0;
    indexToKeepNAS = manual_label ~= 0;
    TP = sum(predicted_label(indexToKeepAS) == 0);                              % wach in T und P
    TN = sum(predicted_label(indexToKeepNAS) ~= 0);
    FP = sum(predicted_label(indexToKeepNAS) == 0);                             % nicht wach in T aber W in P
    FN = sum(predicted_label(indexToKeepAS) ~= 0);
    ACCAS = (TP+TN)/(TP+TN+FP+FN);
    SensAS = TP/(TP+FN);
    SpesAS = TN/(TN+FP);
    

    indexToKeepQS = manual_label == -1;
    indexToKeepNQS = manual_label ~= -1;
    TP = sum(predicted_label(indexToKeepQS) == -1);                             % wach in T und P
    TN = sum(predicted_label(indexToKeepNQS) ~= -1);
    FP = sum(predicted_label(indexToKeepNQS) == -1);                            % nicht wach in T aber W in P
    FN = sum(predicted_label(indexToKeepQS) ~= -1);
    ACCQS = (TP+TN)/(TP+TN+FP+FN);
    SensQS = TP/(TP+FN);
    SpesQS = TN/(TN+FP);
    

    name{end+1} = fileList(fileIdx).name(1:end-12);
    name = strrep(name, '_', ''); % Entfernt "_" im Patientennamen

    Acc_Pat = [Acc_Pat, (ACCW+ACCAS+ACCQS)/3;];
    Sens_Pat = [Sens_Pat , (SensW+SensAS+SensQS)/3];
    Spez_Pat = [Spez_Pat , (SpesW+SpesAS+SpesQS)/3;];
    Kappa_Pat = [Kappa_Pat , kappa_p];
end

evaluation = table(name', Acc_Pat', Sens_Pat', Spez_Pat', Kappa_Pat',   'VariableNames', {'Name','Accuracy', 'Sensitivity', 'Specificity', 'Kappa'});
