%% Extract piezo feature
load 'D:\Sleep_classification\2_processedData\preprocessed_data.mat', 'daten' ;     % Input: preprocessed data
output = 'D:\Sleep_classification\';                                                % Output: piezo feature

%%
fieldNames = fieldnames(daten);
samplingrate = 200;

for i = 1:length(fieldNames) 
    if contains(fieldNames{i}, 'Piezo')
        rawSignal = daten.(fieldNames{i});
        ecgSignal = daten.(fieldNames{i-1});
        hypnogram1 = daten.(fieldNames{i+1});
        hypnogram2 = daten.(fieldNames{i+2});
        hypnogram3 = daten.(fieldNames{i+3});
        plot_title = fieldNames(i);     % "Piezo";

        % time vectors
        sr_Piezo = 200;
        sr_SP = 1/30;
        t_Piezo = 0:1/sr_Piezo:(length(rawSignal)-1)/sr_Piezo;       
        t_SP = 0:1/sr_SP:(length(hypnogram1)-1)/sr_SP;

        % Movement Detection Piezo
        % Spectrogram Piezo
        fs = sr_Piezo;
        win_overlap = 0.5;                                      % overlap in %
        freqs = round(logspace(log10(0.25),log10(100),50),1);   % 1-100 Hz
        freqs = unique(freqs);                                  % spectrogram doesn't accept double values (doppelte Werte)

        % define Paramaters found in ParameterEvaluierung
        nSec = 1;               % length of sliding window in s for spectrogram
        bands_move = [2 25];    % frequency range for bandpower calculation
        std_factor = 4.75;      % factor for calculating the threshold for movement detection
        time_window_s = 5;      % size of moving window in sec, to get rid of gaps in movement detection

        % same code as above, run again with manually set parameters
        Spec_fs = 1/(nSec-nSec*win_overlap); % Time resolution of Spectrogram in Hz
        [s,freqs,t] = spectrogram(rawSignal,nSec*fs,nSec*fs*win_overlap,freqs,fs);
        s = abs(s); % da s auch imaginäre Zahlen enthält
        s_log = 20*log10(abs(s)); % conversion to dB as RMS-spectrum (root-mean-square spectrum)

        % Bandpower
        bandpower_Piezo.bp = bandpower(s, freqs, bands_move, 'psd');    % absolute bandpower across time
        bandpower_Piezo.baseline = prctile(bandpower_Piezo.bp, 25);     % 25th percentile as constant baseline bandpower
        bandpower_Piezo.mad = mad(bandpower_Piezo.bp, 1);               % median absolute deviation (Median ist weniger anfällig bezüglich Ausreißer)
        bandpower_Piezo.threshold = bandpower_Piezo.baseline+std_factor*bandpower_Piezo.mad;

      
        % Movement
        movement = bandpower_Piezo.bp > bandpower_Piezo.threshold;  % bandpower > threshold
        movement = movmax(movement, 5);                             % extend each 1 by 2 datapoints in each direction
        time_window_dp = Spec_fs*time_window_s;                     % time window in data points (dp)
        win_movsum = 2*floor(time_window_dp/2)+1;                   % math to get an odd number for the moving windows
        half_time_window = (win_movsum - 1) / 2;                    % half of it
        m = movsum(movement, win_movsum);                           % moving sum calculation
        movement_smoothed = m >= half_time_window;                    % decision
        clear m

        sr_mov = 1/0.5;
        t_Movement = 0:1/sr_mov:(length(movement_smoothed)-1)/sr_mov;
        % t_vec_Spectrogram_dur = duration((seconds(nSec/2):seconds(1/Spec_fs):length(rawSignal)-seconds(1/Spec_fs)), 'Format','hh:mm:ss'); % 'Format','hh:mm:ss.SSS' to show ms as well
        mov_factor = 1;
        
        window_movement_long = 1200; %600;%1200;
        window_movement_short = 120;
        
        movementInWindow = movsum(movement_smoothed,window_movement_long);

        x = movementInWindow';
        % histogram(bandpower_EEG_smoothed, x)
        [N,edges] = histcounts(x);
        pd = fitdist(x, 'kernel');
        x_values = min(x):0.1:max(x);
        y = pdf(pd, x_values);
        x_localmax = islocalmax(y, 'MaxNumExtrema', 2);
        x_localmin = islocalmin(y, 'MaxNumExtrema', 2);
        y_max = y(x_localmax);
        y_min = y(x_localmin);
        x_max = x_values(x_localmax);
        x_min = x_values(x_localmin);
        movement_baseline = prctile(window_movement_short, 25);
        movement_mad = mad(window_movement_short, 1); % bew_hypnomedian absolute deviation (median is less susceptible to outliers)
        movement_threshold1 = movement_baseline+mov_factor*movement_mad;
        movement_threshold2 = movement_baseline+mov_factor/3*movement_mad;
        movement_threshold1 = 1;
        movement_threshold2 = 1;
        movementInWindow2 = movsum(movement_smoothed,window_movement_short); 
        movement_baseline2 = prctile(movementInWindow2, 25);
        movement_mad2 = mad(movementInWindow2, 1); 
        movement_threshold12 = movement_baseline2+mov_factor*movement_mad;
        movement_threshold22 = movement_baseline2+mov_factor/3*movement_mad;


        movement_hypno = NaN(size(movementInWindow));
        for k=1:length(movementInWindow)
            if movementInWindow(k)> movement_threshold1
                movement_hypno(k) = 1;
            elseif movementInWindow(k) > movement_threshold2
                movement_hypno(k) = 0;
            else
                movement_hypno(k) = -1;
            end
        end
          movement_hypno2 = NaN(size(movementInWindow2));
        for k=1:length(movementInWindow2)
            if movementInWindow2(k)> movement_threshold12
                movement_hypno2(k) = 1;
            elseif movementInWindow2(k) > movement_threshold22
                movement_hypno2(k) = 0;
            else
                movement_hypno2(k) = -1;
            end
        end

        % combi hypno
        movement_hypno_comb = NaN(size(movementInWindow));
        for k=1:length(movementInWindow)
            if movement_hypno(k)==1 & movement_hypno2(k)==0
                movement_hypno_comb(k) = 1;
            elseif movement_hypno(k)==1 & movement_hypno2(k)==-1
                movement_hypno_comb(k) = 0;
            elseif movement_hypno(k)== 0 & movement_hypno2(k)==-1
                movement_hypno_comb(k) = 0;
            elseif movement_hypno(k)== 0 & movement_hypno2(k)==0
                movement_hypno_comb(k) = 1;
            elseif movement_hypno(k)== -1 & movement_hypno2(k)==0
                movement_hypno_comb(k) = 0;    
            elseif movement_hypno(k)== -1 & movement_hypno2(k)==-1
                movement_hypno_comb(k) = -1;    
            else
                movement_hypno_comb(k) = 0.8;
            end
        end

        % average time between peaks
        original_movement = movement;
        new_movement = inf(size(original_movement));

        % Go through the array from left to right
        dist = inf;
        for k = 1:length(original_movement)
            if original_movement(k) == 1
                dist = 0;
            end
            new_movement(k) = dist;
            dist = dist + 1;
        end

        % Move from right to left through the array
        dist = inf;
        for k = length(original_movement):-1:1
            if original_movement(k) == 1
                dist = 0;
            end
            new_movement(k) = min(new_movement(k), dist);
            dist = dist + 1;
        end


        % max dist
        original_array = new_movement;

        % array length
        n = length(original_array);
        transformed_array = zeros(1, n);
        k = 1;
        while k <= n
            if original_array(k) == 0
                transformed_array(k) = 0;
                k = k + 1;
            else
                % Find the end of area without zero
                start_idx = k;
                while k <= n && original_array(k) ~= 0
                    k = k + 1;
                end
                end_idx = k - 1;

                % Determine the maximum value in this area
                max_val = max(original_array(start_idx:end_idx));

                for j = start_idx:end_idx
                    transformed_array(j) = max_val;
                end
            end
        end


        % dists to the next movement
        dist_movement = NaN(size(movementInWindow));
        j = 0;
        for k = 0:length(movementInWindow)-1
            if movement(k+1)== 1
                dist_movement(k+1) = 0;
                j = 0;
            else                 
                j = j+1;                
                dist_movement(k+1) = j;
            end            
        end

        % Find the indices of the ones in the array
        indices_ones1 = find(movement == 1);

        % Create a new array one position larger
        indices_ones = zeros(1, length(indices_ones1) + 1);

        % Insert the 0 in the first position
        indices_ones(2:end) = indices_ones1;
        % Initialisation of the array for the dists with zeros
        dist_array = zeros(size(movement));

        % Calculation of the dists between the ones and storage in dist_array
        for k = 2:length(indices_ones)
            dist_array(indices_ones(k)) = indices_ones(k) - indices_ones(k-1) - 1;
        end

        dist_time = movsum(dist_movement,window_movement_long);
        dist_time2 = movsum(dist_movement,window_movement_short);
        distance_time = movsum(dist_array,window_movement_long);
        distance_time2 = movsum(dist_array,window_movement_short);
        spaceing = movsum(new_movement,window_movement_long);
        spaceing2 = movsum(new_movement,window_movement_short);
        dist_max = movsum(transformed_array,window_movement_long);
        dist_max2 = movsum(transformed_array,window_movement_short);
        
        % Create feature table 30s    
        Label1 = daten.(fieldNames{i+1});
        Label2 = daten.(fieldNames{i+2});
        Label3 = daten.(fieldNames{i+3});
        Label4 = daten.(fieldNames{i+4});
       
        original_array = bandpower_Piezo.bp;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            % z = z*2+1;
            updated_array = [updated_array, z];
        end
        Bandpower = updated_array';
        PiezoFeature_Table = table(Bandpower);

        original_array = movement_smoothed;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Movement = updated_array';

        original_array = movementInWindow;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Activity_l = updated_array';

        original_array = movementInWindow2;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Activity_s = updated_array';

        original_array = movement_hypno;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.hypno1 = updated_array';

        original_array = movement_hypno2;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.hypno2 = updated_array';
        
        original_array = movement_hypno_comb;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.hypno_comb = updated_array';

        original_array = dist_movement;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Inactivity_l = updated_array';

        original_array = dist_array;
        updated_array = [];
        for x = 30: 60: length(original_array)-35   % Value is extracted every 60 entries --> every 30s
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Inactivity_2 = updated_array';

        original_array = new_movement;
        updated_array = [];
        for x = 30: 60: length(original_array)-35   
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Inactivity_3 = updated_array';

        original_array = dist_time;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Inactivity_1_time = updated_array';

        original_array = dist_time2;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Inactivity_1_time2 = updated_array';

        original_array = distance_time;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Inactivity_2_time = updated_array';

        original_array = distance_time2;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Inactivity_2_time2 = updated_array';

        original_array = spaceing;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Inactivity_3_time = updated_array';

        original_array = spaceing2;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.Inactivity_3_time2 = updated_array';

        original_array = dist_max;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.maxstime1 = updated_array';

        original_array = dist_max2;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.maxtime2 = updated_array';

        original_array = transformed_array;
        updated_array = [];
        for x = 30: 60: length(original_array)-35
            z = original_array(x);
            updated_array = [updated_array, z];
        end
        PiezoFeature_Table.dist_max = updated_array';

        while length(updated_array')< length(Label1)
            Label1 = Label1(1:end-1);
            Label2 = Label2(1:end-1);
            Label3 = Label3(1:end-1);
            Label4 = Label4(1:end-1);
        end
        PiezoFeature_Table.Sleepstage1 = Label1;
        PiezoFeature_Table.Sleepstage2 = Label2;
        PiezoFeature_Table.Sleepstage3 = Label3;
        PiezoFeature_Table.Sleepstage4 = Label4;

        folderPath = output;

        % Assuming fieldNames{i} is a valid file name without the extension
        fileName = [fieldNames{i},'_feature', '.csv'];

        % Create the full file path
        fullPath = fullfile(folderPath, fileName);
        writetable(PiezoFeature_Table, fullPath);
    end
end




%% Plot
tiledlayout(5,1)

ax1 = nexttile;
plot(t_Piezo ,rawSignal);
ylabel('Piezo');
title(plot_title,'Interpreter', 'none');
ylim([-1.1 1.1]);

%ECG information
% ax2 = nexttile;
% plot(t_Piezo ,correctes);
% ylabel('filtered Piezo');
% title(plot_title,'Interpreter', 'none');
% ylim([-1.1 1.1]);

% ax3 = nexttile;
% plot(t_Piezo, ecgSignal)

ax4 = nexttile;
hold on
p1 = plot(t_Movement, movementInWindow);
p2 = plot(t_Movement, movementInWindow2);
yline(movement_threshold1);
yline(movement_threshold2);
hold off
glaettung_1 = [num2str(window_movement_long/2/60),' min'];
glaettung_2 = [num2str(window_movement_short/2/60),' min'];
legend([p1 p2],{glaettung_1,glaettung_2}, Location="northeast")
legend('boxoff')
ylabel('Aktivität'); % Anteil der Bewegung in festgelegtem Zeitfenster

ax5 = nexttile;
stairs(t_SP,hypnogram2',Color=[0.4660 0.6740 0.1880], LineWidth=0.9);
ylabel({'Hypnogramm';' manuell'});
ylim([-1.1 1.1]);yticks([-1, 0, 1]); yticklabels({'QS', 'AS', 'W'});

xlabel('Zeit (s)');

linkaxes(findall(gcf, 'type', 'axes'), 'x');


%% PLot: Piezo Feature
tiledlayout(7,1,  'TileSpacing', 'compact')

% Plot 1: ECG signal
ax1 = nexttile;
plot(t_Piezo, rawSignal);
ylabel('Piezo [V]', 'FontSize', 15);
set(gca, 'FontSize', 15, 'XTick', []); % remove ticks

% Plot 2: Spectrogram
ax2 = nexttile;
imagesc(t_Piezo, freqs, 10*log10(abs(s)));
set(gca, 'YTick', freqs(1:10:end), 'CLim', [-60 40], 'FontSize', 16, 'XTick', []); % Colorbar
axis xy;
colorbar;
ylabel('Frequency', 'FontSize', 16);

yticks([0.3, 1.6, 5.3, 18.1, 61.3]);
yticklabels({' ', '','' , '18.1', '61.3'});

movement_smoothed_plot=double(movement_smoothed);
movement_smoothed_plot(movement_smoothed_plot ==0) = NaN;
movement_smoothed_plot = movement_smoothed_plot + 300;
% Plot 3: Bandpower
ax3 = nexttile;
hold on;
plot(t_Movement, bandpower_Piezo.bp);
plot(t_Movement, movement_smoothed_plot, 'LineWidth', 15, 'Color', [0.8500 0.3250 0.0980 0.5]); 
yline(bandpower_Piezo.baseline);
yline(bandpower_Piezo.threshold);
hold off;
ylim([0 340])
ylabel('Bandpower', 'FontSize', 15);
set(gca, 'FontSize', 15, 'XTick', []);

% Plot 4: Movement
% ax4 = nexttile;
% stairs(t_Movement, movement_smoothed);
% ylim([-0.1 1.1]);
% ylabel('Movement', 'FontSize', 15);
% set(gca, 'FontSize', 15, 'XTick', []); 

% Plot 5: Activity
ax4 = nexttile;
hold on;
p1 = plot(t_Movement, movementInWindow);
p2 = plot(t_Movement, movementInWindow2);
yline(movement_threshold1);
yline(movement_threshold2);
hold off;
glaettung_1 = [num2str(window_movement_long/2/60), ' min'];
glaettung_2 = [num2str(window_movement_short/2/60), ' min'];
legend([p1 p2], {glaettung_1, glaettung_2}, 'Location', 'northeast', 'FontSize', 15);
legend('boxoff');
ylabel('Activity', 'FontSize', 15);
set(gca, 'FontSize', 15, 'XTick', []); 

% Plot 6: dist between Movements
ax5 = nexttile;
hold on;
plot(t_Movement, transformed_array);
hold off;
ylabel('dist', 'FontSize', 15);
set(gca, 'FontSize', 15, 'XTick', []); 

% Plot 7: Inactivity over time
ax6 = nexttile;
hold on;
p5 = plot(t_Movement, dist_max, 'Color', [0.6350 0.0780 0.1840]);
p6 = plot(t_Movement, dist_max2);
hold off;
legend([p5 p6], {'10 min', '1 min'}, 'Location', 'northeast', 'FontSize', 15);
legend('boxoff');
ylabel('Inactivity', 'FontSize', 15);
set(gca, 'FontSize', 15, 'XTick', []); 

% Plot 8: Hypnogram
ax7 = nexttile;
stairs(t_SP, hypnogram2', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 0.9);
ylabel('Hypnogram', 'FontSize', 15);
ylim([-1.1 1.1]);
yticks([-1, 0, 1]);
xticks([0, 5000, 10000, 15000]);
xticklabels({'0', '1.5', '3', '4.5'});
yticklabels({'QS', 'AS', 'W'});
xlabel('time (hr)', 'FontSize', 15);
set(gca, 'FontSize', 15); 

% link
linkaxes(findall(gcf, 'type', 'axes'), 'x');



