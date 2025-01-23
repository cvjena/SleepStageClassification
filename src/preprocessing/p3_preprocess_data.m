%% Preprocessing
clear
load('D:\Data_ECG_Piezo.mat', 'daten');              % load concatenated data (ECG, Piezo, Label)
outputfile = 'D:\raw_ECG.fig';                       % output Plot raw signals
outputfile2 = 'D:\processed_ECG.fig';                % output Plot preprocessed signals
output_data = 'D:\preprocessed_data.mat';            % output preprocessed data

%% Plot ECG raw data
originalDaten = daten;
fieldNames = fieldnames(daten);
numECG = sum(contains(fieldNames, 'ECG'));
numSubplotsPerRow = 1;
numRows = ceil(length(fieldNames) / numSubplotsPerRow);

figure;
sr = 200;                                           % sampling rate
sgtitle('ECG raw signal')
tiledlayout(ceil(numECG/3),3)
for i = 1:length(fieldNames)
    if contains(fieldNames{i}, 'ECG')               % field contains "ECG"?
        t_all = 0:1/sr:(length(daten.(fieldNames{i}))-1)/sr;
        ax(i) = nexttile;
        plot(ax(i),t_all,daten.(fieldNames{i}))
        ylabel(ax(i),[fieldNames{i}], 'Interpreter', 'none');
        ylim(ax(i),[-2300,2300]);
    end
end
xlabel('time in s');
linkaxes(findall(gcf, 'Type', 'axes'), 'x');        % link axes

% save
savefig(outputfile);
%% preprocessing

% cut data - manuell bestimmt wenn keine Peaks zu sehen (kein EKG Signal
% vorhanden, bis zu diesem Zeitpunkt wird alles verworfen)
daten.ECG_22 = originalDaten.ECG_22(840000:end); % # epochs * 30 s * 200Hz
daten.Piezo_22 = originalDaten.Piezo_22(840000:end);
daten.Sleepstage_1_22(1:140)=[]; % 140 Epochen werden gelöscht, entspricht 840000 Datenpunkte im ECG
daten.Sleepstage_2_22(1:140)=[];
daten.Sleepstage_3_22(1:140)=[];
daten.Sleepstage_4_22(1:140)=[];
daten.Sleepstage_5_22(1:140)=[];
daten.Sleepstage_6_22(1:140)=[];
daten.Sleepstage_7_22(1:140)=[];

daten.ECG_3 = originalDaten.ECG_3(150000:end);
daten.Piezo_3 = originalDaten.Piezo_3(150000:end);
daten.Sleepstage_1_3(1:25)=[];
daten.Sleepstage_2_3(1:25)=[];
daten.Sleepstage_3_3(1:25)=[];
daten.Sleepstage_4_3(1:25)=[];
daten.Sleepstage_5_3(1:25)=[];
daten.Sleepstage_6_3(1:25)=[];
daten.Sleepstage_7_3(1:25)=[];

daten.ECG_29 = originalDaten.ECG_29(1320000:end);
daten.Piezo_29 = originalDaten.Piezo_29(1320000:end);
daten.Sleepstage_1_29(1:220)=[];
daten.Sleepstage_2_29(1:220)=[];
daten.Sleepstage_3_29(1:220)=[];
daten.Sleepstage_4_29(1:220)=[];
daten.Sleepstage_5_29(1:220)=[];
daten.Sleepstage_6_29(1:220)=[];
daten.Sleepstage_7_29(1:220)=[];


% mirror data (R peaks up) - manuell bestimmt wenn maximaler peak nach oben zeigt
daten.ECG_1  = -originalDaten.ECG_1;
daten.ECG_5  = -originalDaten.ECG_5;
daten.ECG_6  = -originalDaten.ECG_6;
daten.ECG_10 = -originalDaten.ECG_10;
daten.ECG_12 = -originalDaten.ECG_12;
daten.ECG_14 = -originalDaten.ECG_14;
daten.ECG_15 = -originalDaten.ECG_15;
daten.ECG_16 = -originalDaten.ECG_16;
daten.ECG_17 = -originalDaten.ECG_17;
daten.ECG_18 = -originalDaten.ECG_18;
daten.ECG_21 = -originalDaten.ECG_21;
daten.ECG_22 = -daten.ECG_22;           % "daten.ECG_22" weil diese schon gekürzt wurden
daten.ECG_23 = -originalDaten.ECG_23;
daten.ECG_29 = -daten.ECG_29;           % daten.ECG_29 weil diese schon gekürzt wurden

% Plot
figure;
sgtitle('ECG preprocessed')
tiledlayout(ceil(numECG/3),3)
for i = 1:length(fieldNames)
    if contains(fieldNames{i}, 'ECG')

        t_all = 0:1/sr:(length(daten.(fieldNames{i}))-1)/sr;
        ax(i) = nexttile;
        plot(ax(i),t_all,daten.(fieldNames{i}))
        ylabel(ax(i),[fieldNames{i}], 'Interpreter', 'none');
        ylim(ax(i),[-2300,2300]);
        set(gca, 'Position', get(gca, 'Position') + [0, -0.05, 0, 0]);
    end
end
xlabel('time in s');
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

save(output_data, 'daten');
savefig(outputfile2);







