%% Extract ECG fatures with optimizes r peak detection and without HRV tool
clear
load('D:\preprocessed_data.mat', 'daten');     % Input: preprocessed data
output = 'D:\ECG2_feature\';                                 % Output: ECG2 features
output_ECG2_feature = 'D:\ECG2_feature\';

%%
fieldNames = fieldnames(daten);
for ii = 1:3%length(fieldNames) 
    if contains(fieldNames{ii}, 'ECG')      % contains ECG?
        signal = daten.(fieldNames{ii});
        hypno = daten.(fieldNames{ii+2});       % Dataset A
        hypno2 = daten.(fieldNames{ii+3});      % Dataset B
        plot_title = "ECG ";

        % time vectors
        sr_ECG = 200;
        t_ECG = 0:1/sr_ECG:(length(signal)-1)/sr_ECG;
        sr_SP = 1/30;
        t_SP = 0:1/sr_SP:(length(hypno)-1)/sr_SP;

        % fft
        fresult = fft(signal);
        fresult(1:round(length(fresult)*5/sr_ECG))=0;
        fresult(end - round(length(fresult)*5/sr_ECG):end)=0;
        correctes=real(ifft(fresult));

        % peak detection
        [pks,locs] = findpeaks(correctes,t_ECG, 'MinPeakProminence',60,'MinPeakDistance',0.2, 'MaxPeakWidth',0.1703, 'WidthReference','halfprom', 'MinPeakHeight',120);
        pks(end) =[];           
        locs(end)=[];           
        peak_distance = diff(locs);                 % distance between peaks

        
        % remove false positive peaks
        for i = 2:(length(peak_distance)-1)
            try
                if peak_distance(i) <= (peak_distance(i-1)/1.5) && peak_distance(i-1)<=1 && peak_distance(i) <=(peak_distance(i+3)/1.5)
                    % If peak distances are unusually short, the
                    % associated peak and its position are deleted
                    pks(i+1)=[];
                    locs(i+1)=[];
                    peak_distance = diff(locs); 
                    i = i-1;
                    continue; 
                end
            catch ME
                break;
            end
        end

        % calculate spectrogram         
        nSec = 3;
        fs = 200;
        win_overlap = 0.5;
        Spec_fs = 1/(nSec-nSec*win_overlap);                % Time resolution of Spectrogram in Hz
        freqs = round(logspace(log10(1),log10(30),30),1);   % 1-50 Hz
        freqs = unique(freqs);                              % spectrogram doesn't accept double values
        [s,freqs,t_ECG_SG] = spectrogram(signal,nSec*fs,nSec*fs*win_overlap,freqs,fs);
        s = abs(s);
        s_log = 20*log10(abs(s));                           % conversion to dB as RMS-spectrum (root-mean-square spectrum)

        % calculate bandpower
        bands_ECG = [1 30];
        bandpower_ECG = bandpower(s, freqs, bands_ECG, 'psd');
        bandpower_ECG_smoothed = movmedian(bandpower_ECG,5);
        bandpower_ECG_baseline = prctile(bandpower_ECG, 10); % 25th percentile as constant baseline bandpower
        bandpower_ECG_mad = mad(bandpower_ECG, 1);           % median absolute deviation
        ECG_thr_factor = 15;
        bandpower_ECG_threshold = bandpower_ECG_baseline + ECG_thr_factor * bandpower_ECG_mad; % EEG threshold

        % timevector BP
        sr_BP = 1/1.5;
        t_BP = 0:1/sr_BP:(length(bandpower_ECG)-1)/sr_BP;

        % calculate noise
        noise = bandpower_ECG >= bandpower_ECG_threshold;
        noise_plot = zeros(size(noise));
        noise_plot = noise_plot';

        for i = 1:length(noise)
            if noise(i) == 1
                noise_plot(i-3:i+3) = 1000000;
            else
            end
        end

        % searches for points (seconds) at which artefacts occur in order to set these later
        % to NaN
        positions1 = find(noise_plot > 0);
        positions = flip(positions1);
        neu = positions1;
        for i = 2:length(positions)-1
            if positions(i) == positions(i-1)-1 && positions(i) == positions(i+1)+1
                neu(length(positions)-i+1) = [];
            end
        end
        neu = neu*1.5;          % time in s
        neu(2:end+1) = neu;     % Move one position to the right
        neu(1) = 0;             % Insert start position

        % Set peak distances in regions with noise to NaN
        for i = 2:2:length(neu)
            for j = 1:length(locs)-1
                if locs(j) > neu(i) && locs(j) < neu(i+1)
                    peak_distance(j) = NaN;
                end
            end
        end

        % calculate features
        mean_R_dist_l = movmean(peak_distance, 20,'omitnan');   % beachtet NaN's -> EKG Bereiche die aufgrund von Artefakten entfernt wurden
        mean_R_dist_h = movmean(peak_distance, 100,'omitnan');

        median_R_dist_l = movmedian(peak_distance, 20,'omitnan');
        median_R_dist_h = movmedian(peak_distance, 100,'omitnan');

        std_R_dist_l = movstd(peak_distance, 20,'omitnan');
        std_R_dist_h = movstd(peak_distance, 100,'omitnan');

        difference_mean_median = mean(abs(mean_R_dist_l- median_R_dist_h),10,'omitnan');

        coefficent_of_var_l = std_R_dist_l./mean_R_dist_l;
        coefficent_of_var_h = std_R_dist_h./mean_R_dist_h;

        % Set feature in noise areas to NaN
        nan_indices = isnan(peak_distance);
        numOnes = sum(nan_indices == 1, 'all');

        % numNaN = sum(isnan(mean_R_dist_h), 'all');

        mean_R_dist_l(nan_indices) = NaN;
        mean_R_dist_h(nan_indices) = NaN;

        median_R_dist_l(nan_indices) = NaN;
        median_R_dist_h(nan_indices) = NaN;

        std_R_dist_l(nan_indices) = NaN;
        std_R_dist_h(nan_indices) = NaN;

        difference_mean_median(nan_indices) = NaN;

        coefficent_of_var_l(nan_indices) = NaN;
        coefficent_of_var_h(nan_indices) = NaN;

        % use 30s intervals

        imean_L = intervalt(mean_R_dist_l, locs);           %intervalt: lokale Funktion s.u.
        imean_R_dist_h = intervalt(mean_R_dist_h, locs);
        imedian_R_dist_l = intervalt(median_R_dist_l, locs);
        imedian_R_dist_h = intervalt(median_R_dist_h, locs);
        istd_R_dist_h = intervalt(std_R_dist_h, locs);
        istd_R_dist_l = intervalt(std_R_dist_l, locs);
        idiff = intervalt(difference_mean_median, locs);
        icoff_l = intervalt(coefficent_of_var_l, locs);
        icoff_h = intervalt(coefficent_of_var_h, locs);

        % save feature values

        if length(hypno)> length(imean_L)
            Feature_table = table(hypno(1:end-1), 'VariableNames', {'Sleepstage1'});
        else
            Feature_table = table(hypno, 'VariableNames', {'Sleepstage1'});
        end
        Feature_table.mean_l = imean_L';
        Feature_table.mean_h = imean_R_dist_h';
        Feature_table.median_l = imedian_R_dist_l';
        Feature_table.median_h = imedian_R_dist_h';
        Feature_table.std_h = istd_R_dist_h';
        Feature_table.std_l = istd_R_dist_l';
        Feature_table.diff = idiff';
        Feature_table.cff_l = icoff_l';
        Feature_table.cff_h = icoff_h';

        % save
        folderPath = output;

        fileName = [fieldNames{ii},'_feature', '.csv'];

        % Create the full file path
        fullPath = fullfile(folderPath, fileName);
        writetable(Feature_table, fullPath);

        % Plot 
        figure
        tiledlayout(6,1)

        ax1 = nexttile;
        hold on
        plot(t_ECG,correctes);
        % plot(locs1, pks1, 'ro', 'MarkerSize', 10, 'MarkerEdgeColor','b');
        plot(locs, pks, 'ro', 'MarkerSize', 10);  
        hold off
        ylabel('ECG');

        ax2 = nexttile;
        imagesc(t_ECG, freqs, 10*log10(abs(s)));
        axis xy;
        colorbar;
        ylabel('Spectrogram');
     
        ax3 = nexttile;
        hold on
        plot(t_BP,bandpower_ECG);
        yline(bandpower_ECG_baseline, 'g');
        yline(bandpower_ECG_mad, 'b');
        yline(bandpower_ECG_threshold, 'r');
        %stairs(t_BP, noise_plot, 'b');
        hold off
        ylabel('Bandpower')

        ax4 = nexttile;
        hold on
        % scatter(peak_times(2:end), peak_distances);       % original
        % R-Peak distances
        scatter(locs(2:end), peak_distance);                % reduced distances (NaN)

        % plot(time_array, mean_distance_array, '-o');
        plot(locs(2:end), mean_R_dist_l);
        plot(locs(2:end), mean_R_dist_h);
        plot(locs(2:end), median_R_dist_l);
        plot(locs(2:end), median_R_dist_h);
        legend('distance', 'mean (short)','mean (long)', 'median (short)', 'median (long)');
        hold off
        ylabel('R-Peak distance [s]');
                
        ax5 = nexttile;
        hold on
        plot(locs(2:end), std_R_dist_l);
        plot(locs(2:end), std_R_dist_h);
        plot(locs(2:end), difference_mean_median);
        plot(locs(2:end), coefficent_of_var_l);
        plot(locs(2:end), coefficent_of_var_h);
        legend('std (short)', 'std (long)', '|mean-median|', 'CoV (short)', 'CoV (long)');
        hold off
        ylabel('Std, abs, cof of var');

        ax6 = nexttile;
        hold on
        stairs(t_SP,hypno2);
        hold off
        ylabel('Hypnogram (B)');
        ylim([-1.1 1.1]);yticks([-1, 0, 1]); yticklabels({'QS', 'AS', 'W'});xlabel('Zeit in s');

        linkaxes(findall(gcf, 'type', 'axes'), 'x');
        
        filename = [output_ECG2_feature,(fieldNames{ii}),'_ECG2_feature.fig']; % ECG_1 steht für Patient 1
        savefig(filename);
        clf
    end
end

%% !
% Minimise distances
tiledlayout(5,1, 'TileSpacing', 'tight', 'Padding', 'tight');
fontSize = 18;

% Plot 1
ax1 = nexttile;
hold on
plot(t_ECG, correctes);
% plot(locs1, pks1, 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'b');
plot(locs, pks, 'ro', 'MarkerSize', 10);  
hold off
ylabel('ECG [μV]', 'FontSize', fontSize);
ax1.XColor = 'none'; 
ax1.YAxis.FontSize = fontSize; 
ax1.XAxis.FontSize = fontSize; 

% Plot 2
ax2 = nexttile;
imagesc(t_ECG, freqs, 10*log10(abs(s)));
axis xy;
colorbar;
ylabel('Frequency', 'FontSize', fontSize);
ax2.XColor = 'none'; 
ax2.YAxis.FontSize = fontSize; 
ax2.XAxis.FontSize = fontSize; 

% Plot 3
ax3 = nexttile;
hold on
plot(t_BP, bandpower_ECG);
yline(bandpower_ECG_baseline, 'g');
yline(bandpower_ECG_mad, 'b');
yline(bandpower_ECG_threshold, 'r');
hold off
ylabel('Bandpower', 'FontSize', fontSize);
ax3.XColor = 'none'; 
ax3.YAxis.FontSize = fontSize; 
ax3.XAxis.FontSize = fontSize; 

% Plot 4
ax4 = nexttile;
hold on
scatter(locs(2:end), peak_distance);               
plot(locs(2:end), mean_R_dist_l);
plot(locs(2:end), mean_R_dist_h);
plot(locs(2:end), median_R_dist_l);
plot(locs(2:end), median_R_dist_h);
legend('distance', 'mean (short)', 'mean (long)', 'median (short)', 'median (long)', 'FontSize', fontSize);
legend("boxoff")
hold off
ylabel('RR interval [s]', 'FontSize', fontSize);
ax4.XColor = 'none'; 
ax4.YAxis.FontSize = fontSize; 
ax4.XAxis.FontSize = fontSize; 

% Plot 5
ax5 = nexttile;
hold on
plot(locs(2:end), std_R_dist_l);
plot(locs(2:end), std_R_dist_h);
plot(locs(2:end), difference_mean_median);
plot(locs(2:end), coefficent_of_var_l);
plot(locs(2:end), coefficent_of_var_h);
legend('std (short)', 'std (long)', '|mean-median|', 'CoV (short)', 'CoV (long)', 'FontSize', fontSize);
legend("boxoff")
hold off
ylabel('Feature', 'FontSize', fontSize);
xlabel('time [s]', 'FontSize', fontSize); 
ax5.YAxis.FontSize = fontSize; 
ax5.XAxis.FontSize = fontSize; 

axesHandles = findall(gcf, 'Type', 'axes');

minDistance = 0.1; 

for i = 1:numel(axesHandles)
    pos = axesHandles(i).OuterPosition;
    pos(1) = pos(1) - minDistance; 
    axesHandles(i).OuterPosition = pos;
end

linkaxes(findall(gcf, 'type', 'axes'), 'x');

%%
function mean_distance_array = intervalt(distances, x_values)
start_time = x_values(1);
end_time = x_values(end);
time_interval = 30; 
time_array = start_time:time_interval:end_time;

mean_distance_array = zeros(size(time_array));

% Find the averaged distances for the 30 second intervals
current_peak = 1;
for i = 1:length(time_array)
    while current_peak < length(x_values) && x_values(current_peak + 1) < time_array(i)
        current_peak = current_peak + 1;
    end
    if current_peak < length(distances)
        mean_distance_array(i) = distances(current_peak);
    else
        mean_distance_array(i) = distances(end);
    end
end
end
