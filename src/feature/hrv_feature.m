%% extract hrv feature from ECG signal
clear
load('D:\preprocessed_data.mat', 'daten'); % Input: preprocessed data
output = 'D:\ECG_feature\';                              % Output: ECG feature
output_peak_detection = 'D:\peak_detection\';
output_peak_distance='D:\peak_distance\';
% add Heart Rate Variability Feature Set
% (https://de.mathworks.com/matlabcentral/fileexchange/109895-heart-rate-variability-feature-set)
%hrv.EAR: ändern: L_overlap = 90; % in percentage damit es alle 30s ein
%feature gibt
% kein Peak künstlich nach 2s ohne Detektion einfügen! --> hrv.EAR: max_rr_interval = 2;
%%
fieldNames = fieldnames(daten);
sr = 200;                                   % sampling rate

for i = 1:length(fieldNames) 
    if contains(fieldNames{i}, 'ECG')       % if field contains "ECG"
        t_ECG = 0:1/sr:(length(daten.(fieldNames{i}))-1)/sr;
        
        % fft
        fresult = fft(daten.(fieldNames{i}));
        fresult(1:round(length(fresult)*5/sr)) = 0;
        fresult(end - round(length(fresult)*5/sr):end) = 0;
        correctes = real(ifft(fresult));

        label1 = daten.(fieldNames{i+2});   % Dataset A
        label2 = daten.(fieldNames{i+3});   % Dataset B
        
        % r peak detection
        [pks,locs] = findpeaks(correctes,t_ECG, 'MinPeakProminence',60,'MinPeakDistance',0.2, 'MaxPeakWidth',0.1703, 'WidthReference','halfprom', 'MinPeakHeight',120);
        pks(end) =[];                   % remove last entry
        locs(end)=[];                   % remove last entry
        peak_distance = diff(locs);     % calculate distance between peaks

        % remove false pos. peaks
        for j = 2:(length(peak_distance)-1)
            try
                if peak_distance(j) <= (peak_distance(j-1)/1.5) && peak_distance(j-1)<=1 && peak_distance(j) <=(peak_distance(j+3)/1.5)
                    % If peak distances are unusually short, the
                    % associated peak and its position are deleted
                    pks(j+1)=[];
                    locs(j+1)=[];
                    peak_distance = diff(locs); % Update distances
                    j = j-1;
                    continue; 
                end
            catch ME
                break;
            end
        end       
        

        time_vector = 1:length(locs);
        
        tiledlayout(3,1)
        ax1 = nexttile;
        plot(t_ECG,correctes);
        hold on;
        plot(locs, pks, 'ro', 'MarkerSize', 10);  
        ylabel('ECG signal')
        ax2 = nexttile;
        scatter(locs(2:end), peak_distance);
        ylabel({'Distance between',' R-Peaks (s)'});
        
        disp(fieldNames{i})
        teile = split(fieldNames{i}, '_');
        letzter_teil = teile{end};
       
        timeNumeric = 0:30:(height(daten.(['Sleepstage_1_' num2str(letzter_teil)]))-1)*30;
        ax3 = nexttile;

        stairs(timeNumeric,daten.(['Sleepstage_2_' num2str(letzter_teil)]))
        ylabel('Hypnogram');
        ylim([-1.1 1.1]);yticks([-1, 0, 1]); yticklabels({'QS', 'AS', 'W'});xlabel('Zeit in s');
        
        linkaxes(findall(gcf, 'type', 'axes'), 'x');

        filename = [output_peak_distance,(fieldNames{i}),'_peakDistance.fig'];
        % Speichern des Plots
        savefig(filename);
        % close(gcf);

        %         %% für plot in bpm
        %         tiledlayout(3,1)
        %         ax1 = nexttile;
        %         plot(t_ECG,correctes);
        %         hold on;
        %         plot(locs, pks, 'ro', 'MarkerSize', 10);  % Peaks markieren
        %         ylabel('ECG signal')
        %         ax = gca; % Get Current Axis
        %
        % % Schriftgröße der Achsenbeschriftungen ändern
        % % ax.XLabel.String = 'X-Achse';
        % % ax.XLabel.FontSize = 20; % Schriftgröße der x-Achsenbeschriftung
        % ax.YLabel.String = 'ECG signal';
        % ax.YLabel.FontSize = 20; % Schriftgröße der y-Achsenbeschriftung
        %
        % % Schriftgröße der Achsenticks ändern
        % ax.XAxis.FontSize = 15; % Schriftgröße der x-Achsenticks
        % ax.YAxis.FontSize = 15; % Schriftgröße der y-Achsenticks
        %
        %         ax2 = nexttile;
        %         % Plot der Abstände über die Zeit
        %         scatter(locs(2:end), 60./peak_distance);
        %         ylabel({'Heart rate [bpm]'});
        %         % Achsenbeschriftungen hinzufügen
        % ax = gca; % Get Current Axis
        %
        % % Schriftgröße der Achsenbeschriftungen ändern
        % % ax.XLabel.String = 'X-Achse';
        % % ax.XLabel.FontSize = 20; % Schriftgröße der x-Achsenbeschriftung
        % ax.YLabel.String = 'Heart rate [bpm]';
        % ax.YLabel.FontSize = 20; % Schriftgröße der y-Achsenbeschriftung
        %
        % % Schriftgröße der Achsenticks ändern
        % ax.XAxis.FontSize = 15; % Schriftgröße der x-Achsenticks
        % ax.YAxis.FontSize = 15; % Schriftgröße der y-Achsenticks
        %
        % title(fieldNames{i}, Interpreter="none");
        % % hypno
        %         disp(i)
        %         disp(fieldNames{i})
        %         teile = split(fieldNames{i}, '_');
        %         letzter_teil = teile{end};
        %         disp(letzter_teil);
        %         timeNumeric = 0:30:(height(daten.(['Sleepstage_1_' num2str(letzter_teil)]))-1)*30;
        %         ax3 = nexttile;
        %
        %
        %         stairs(timeNumeric,daten.(['Sleepstage_2_' num2str(letzter_teil)]))
        %         ylabel('Hypnogram');
        %         ylim([-1.1 1.1]);yticks([-1, 0, 1]); yticklabels({'QS', 'AS', 'W'});xlabel('time [s]');
        %         ax = gca; % Get Current Axis
        %
        % % Schriftgröße der Achsenbeschriftungen ändern
        % ax.XLabel.String = 'time [s]';
        % ax.XLabel.FontSize = 20; % Schriftgröße der x-Achsenbeschriftung
        % ax.YLabel.String = 'Hypnogram';
        % ax.YLabel.FontSize = 20; % Schriftgröße der y-Achsenbeschriftung
        %
        % % Schriftgröße der Achsenticks ändern
        % ax.XAxis.FontSize = 15; % Schriftgröße der x-Achsenticks
        % ax.YAxis.FontSize = 15; % Schriftgröße der y-Achsenticks
        %
        %         linkaxes(findall(gcf, 'type', 'axes'), 'x');

        %
        figure;
        plot(t_ECG,correctes);
        hold on;
        plot(locs, pks, 'ro', 'MarkerSize', 10);  
        
        xlabel('Zeit in s');
        ylabel('Amplitude');
        title(fieldNames{i}, Interpreter="none");
        filename = [output_peak_detection,(fieldNames{i}),'_peakDetection.fig'];
        savefig(filename);
        close(gcf);
        
        % create struct, Input HRV tool
        myStruct.code = 'ID_1';                  % ID patient, in dem Fall immer die gleiche, weil das HRV Tool pro Patient aufgerufen wird
        myStruct.rr_interval = peak_distance;    % rr distances
        myStruct.rr_peaks = locs;                % peak position

        % HRV tool
        [hrv_feats_tb, hrv_feats_epochs_tb] = hrv_features(myStruct); 
       
        string = fieldNames(i);
        String = cell2mat(string);
        parts = strsplit(String, '_');
        lastPart = parts{end};

        tabelle = hrv_feats_epochs_tb;

        for r = 1:4 % 4 leerzeilen einfügen damit feature an epoche 5 =2,5min steht (zentriertes Fenster)
            hrv_feats_epochs_tb = vertcat(array2table(NaN(1, width(hrv_feats_epochs_tb)), 'VariableNames', hrv_feats_epochs_tb.Properties.VariableNames), hrv_feats_epochs_tb);
        end
        
        verbunden = ['Sleepstage_1_', lastPart];
        anzahl_Sleepstage = numel(daten.(verbunden));
       
        anzahl_Zeilen= size(hrv_feats_epochs_tb.mean_NN);
        anzahl_Zeilen = anzahl_Zeilen(1);

        
        differenz = anzahl_Sleepstage - anzahl_Zeilen;
        if anzahl_Sleepstage > anzahl_Zeilen %&& differ <= 5
            for j = 1:differenz
                hrv_feats_epochs_tb = [hrv_feats_epochs_tb; array2table(NaN(1, width(tabelle)), 'VariableNames', tabelle.Properties.VariableNames)];
            end
        end
      
      
        hrv_feats_epochs_tb = removevars(hrv_feats_epochs_tb, 'baby_ID');

        for l = 1:size(hrv_feats_epochs_tb, 2)
            spalteName = hrv_feats_epochs_tb.Properties.VariableNames{l};
            spalte = hrv_feats_epochs_tb.(spalteName);
            werte = daten.(verbunden);
            hrv_feats_epochs_tb.Sleepstage1 = werte;           % Dataset A
        end
        folderPath = output;
        fileName = [fieldNames{i},'_feature', '.csv'];
        fullPath = fullfile(folderPath, fileName);
        hrv_feats_epochs_tb.Sleepstage2 = label2;       % Dataset B
        writetable(hrv_feats_epochs_tb, fullPath);
    end
end
