%-----------------------------Loading ECG Signal---------------------------

%moving to the directory where the dataset is present
cd('D:\Iqra\.DEGREE\SEM 6\DSP\Project\mit-bih-arrhythmia-database-1.0.0'); 

%loading the ecg signal
[ecgSignal, Fs, t] = rdsamp('100',1); % Fs = 360 Hz for MIT-BIH
%as the total signall isof30 min so we just laoded it for 3 sec for better
%understanding
signal_duration = 3; 
sample_extracted = signal_duration * Fs; % Samples = time × sampling freq

%slicing the signal into 3 seconds
ecgSignal = ecgSignal(1:sample_extracted);
t = t(1:sample_extracted);

ecg = ecgSignal(:,1);                  %extracting the one channel from complete signal

%plot clean signall
figure;
plot(t, ecg);
title('Clean ECG Signal');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');

%%
%-------------------------NOISE ADDITION-----------------------------
%Powerline interference (50Hz)
powerline_noise = 0.2*sin(2*pi*50*t);

%Baseline wander (0.3Hz)
baseline_wander = 0.5*sin(2*pi*0.3*t);

%EMG noise (muscle artifact)
emg_noise = 0.15*randn(size(t));

%combining all noises to create a single noisy signal
noisy_ecg = ecg + powerline_noise + baseline_wander + emg_noise;

figure;                    %all the noisy signals seperately
subplot(3,1,1);
plot(t, (powerline_noise + ecg));
title('ECG with Powerline Noise of 50Hz');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');

subplot(3,1,2);
plot(t, (baseline_wander + ecg));
title('ECG with Baseline Wander Noise of 0.3Hz');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');

subplot(3,1,3);
plot(t, (emg_noise + ecg));
title('ECG with EMG Noise');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');

%plotting the clean and noisy signal that combines all the noises to
%compare
figure;
subplot(2,1,1);
plot(t, ecg);
title('Clean ECG');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');

subplot(2,1,2);
plot(t, noisy_ecg);
title('Noisy ECG (50Hz + Baseline + EMG noise)');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');

%%
%-------------------------FREQUENCY ANALYSIS--------------------------------
n = length(ecg);
f = (0:n-1)*(Fs/n);

%computing FFTs for both clean n noisy signal
clean_fft = abs(fft(ecg));
noisy_fft = abs(fft(noisy_ecg));

figure;
subplot(2,1,1);
plot(f, clean_fft);
xlim([0,100]);
title('Clean ECG Spectrum');
xlabel('Frequency (Hz)'); 
ylabel('Magnitude');

subplot(2,1,2);
plot(f, noisy_fft);
xlim([0,100]);
title('Noisy ECG Spectrum');
xlabel('Frequency (Hz)'); 
ylabel('Magnitude');

%Spectrogram-Shows how the frequency content of a signal changes over time.
window_length = 64;     % Smaller window due to short signal
noverlap = 32;
nfft = 128;

figure;
subplot(2,1,1);
spectrogram(ecg, window_length, noverlap, nfft, Fs, 'yaxis');
title('Spectrogram of Clean ECG Signal');
ylim([0 50]); % Focus on important ECG frequencies
colorbar;

subplot(2,1,2);
spectrogram(noisy_ecg, window_length, noverlap, nfft, Fs, 'yaxis');
title('Spectrogram of Noisy ECG Signal');
ylim([0 50]); % Focus on important ECG frequencies
colorbar;

%Short-Time Fourier Transform (STFT)-Same as spectrogram, but sometimes plotted differently (higher customization).
[s_e, f_e, t_e] = stft(ecg, Fs, 'Window', hann(256), 'OverlapLength', 200, 'FFTLength', 512);
[s_n, f_n, t_n] = stft(noisy_ecg, Fs, 'Window', hann(256), 'OverlapLength', 200, 'FFTLength', 512);

% Plot magnitude
figure;
subplot(2,1,1);
imagesc(t_e, f_e, abs(s_e));
axis xy;
title('Time-Frequency Plot (STFT) of Clean ECG');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0 50]);
colorbar;

subplot(2,1,2);
imagesc(t_n, f_n, abs(s_n));
axis xy;
title('Time-Frequency Plot (STFT) of Noisy ECG');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0 50]);
colorbar;


%% 
%-------------------------FILTER DESIGN AND APPLICATION--------------------

%using IIR eliptic filters for better performance
% notch = NF_El(noisy_ecg);           
% highpass = HF_El(notch);               
% filtered_ecg = LF_El(highpass);        

%using IIR eliptic filters as objects for better performance
notch = filter(NF_EL_O,noisy_ecg);           
highpass = filter(HF_EL_O,notch);              
filtered_ecg = filter(LF_EL_O,highpass);        

%HF_Filter, LF_Fillter1, Notch_Filter
%ecg signal after each filter
figure;
subplot(3,1,1);
plot(t, notch);
title('Notch Filter');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');

subplot(3,1,2);
plot(t, highpass);
title('High Pass Filter');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');

subplot(3,1,3);
plot(t, filtered_ecg);
title('Low Pass Filter');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');

%comparison ecg signals
figure;
plot(t, noisy_ecg,'b');
hold on;
plot(t, filtered_ecg,'g');
hold on;
plot(t, ecg,'r');
legend('Noisy ECG', 'Filtered ECG','Clean ECG');
title('Static Filter Comparison');
xlabel('Time (s)'); 
ylabel('Amplitude (mV)');


%%
%--------------------ADAPTIVE FILTER IMPLEMENTATION------------------------
import dsp.*

N = 64;               %length of filter
mu = 0.0001;          %step size

%reference signal for the removal of noise
ref = sin(2*pi*50*t);  
ref = ref(:);                            %column vector
ref = [zeros(N-1, 1); ref];  % Pad the reference signal to match filter order

%checking the length
n = length(noisy_ecg);
if length(ref) > n
    ref = ref(1:n);  
elseif length(ref) < n
    ref = [ref; zeros(n - length(ref), 1)];  
end

%highpass = filter(HP_Filter,noisy_ecg);               
%filtered_ecg_ad = LF_El(noisy_ecg); 

%LMS filter object
lmsFilter = dsp.LMSFilter('Length', N, 'StepSize', mu);
[y_p, e_p] = lmsFilter(ref, noisy_ecg);               %applying the LMS filter to reduce powerline

t = (0:length(e_p)-1)/Fs; 

%baseline reference: low-frequency random drift
rng(1);                                                 %Reproducibility
ref_b = lowpass(randn(size(t)), 0.5, Fs);  
ref_b = ref_b(:);              
e_p = e_p(:);  

[y, e] = lmsFilter(ref_b, e_p);                        %applying lms to reduce vaseline wander

%Median filter to remove small artifacts
adap_filtered_ecg = medfilt1(e, 5);  

figure;
plot(t, noisy_ecg, 'b');
hold on;
plot(t, adap_filtered_ecg, 'g');  
hold on;
plot(t,ecg,'r');
%xlim([0.1,0.15]);
legend('Noisy ECG', 'Filtered ECG','Clean ECG');
xlabel('Time (s)');
ylabel('Amplitude');
title('LMS Adaptive Filtering Using Built-in Function');


%%
%-------------------QUANTITARIVE METRICS--------------------------------
%ensuring all signals are column vectors
clean_ecg = ecg(:);
noisy_ecg = noisy_ecg(:);
filtered_ecg = filtered_ecg(:);
adap_filtered_ecg = adap_filtered_ecg(:);

%finding minimum length to set all signals to it
min_len = min([length(clean_ecg), length(noisy_ecg), length(filtered_ecg), length(adap_filtered_ecg)]);
clean_ecg = clean_ecg(1:min_len);
filtered_ecg = filtered_ecg(1:min_len);
adap_filtered_ecg = adap_filtered_ecg(1:min_len);

%SIGNAL DISTORTION- snr(signal,noise)
initial_snr = snr(clean_ecg,(noisy_ecg - clean_ecg));
static_snr = snr(clean_ecg,(filtered_ecg- clean_ecg));
adaptive_snr = snr(clean_ecg,(adap_filtered_ecg- clean_ecg));

fprintf('-------------------------------------------------------------------\n');
fprintf('SNR Results:\n');
fprintf('Initial SNR: %.2f dB\n', initial_snr);
fprintf('Static Filter SNR: %.2f dB\n', static_snr);
fprintf('Adaptive Filter SNR: %.2f dB\n', adaptive_snr);

%RMSE
rmse_static = sqrt(mean((clean_ecg - filtered_ecg).^2));
rmse_adaptive = sqrt(mean((clean_ecg - adap_filtered_ecg).^2));

fprintf('\n-------------------------------------------------------------------');
fprintf('\nRMSE Results:\n');
fprintf('Static Filter RMSE: %.4f mV\n', rmse_static);
fprintf('Adaptive Filter RMSE: %.4f mV\n', rmse_adaptive);
%%
%R-Peak Detection Comparison
min_peak_distance = round(0.2 * Fs);  %min distance ~200ms (~300 bpm)

[~, loc_clean] = findpeaks(clean_ecg, 'MinPeakProminence', 0.3);
[~, loc_static] = findpeaks(filtered_ecg, 'MinPeakProminence', 0.3);
[~, loc_adaptive] = findpeaks(adap_filtered_ecg, 'MinPeakProminence', 0.3);

%tolerance
tolerance = round(0.1 * Fs);

acc_static = R_peak_Accuracy(loc_clean, loc_static, tolerance);
acc_adaptive = R_peak_Accuracy(loc_clean, loc_adaptive, tolerance);

fprintf('\n-------------------------------------------------------------------');
fprintf('\nR-Peak Detection Accuracy:\n');
fprintf('Static Filter Accuracy: %.2f%%\n', acc_static);
fprintf('Adaptive Filter Accuracy: %.2f%%\n', acc_adaptive);

figure;
plot(clean_ecg); hold on;
plot(loc_clean, clean_ecg(loc_clean), 'ko');  % Reference peaks
plot(loc_static, clean_ecg(loc_static), 'gx');  % Static peaks
plot(loc_adaptive, clean_ecg(loc_adaptive), 'r+');  % Adaptive peaks
legend('ECG', 'Clean', 'Static', 'Adaptive');
title('R-peak Comparison');


function accuracy = R_peak_Accuracy(ref_peaks, test_peaks, tolerance)
    matched = false(size(test_peaks));  % Track matched test peaks
    correct = 0;

    for i = 1:length(ref_peaks)
        diffs = abs(test_peaks - ref_peaks(i));
        [min_diff, idx] = min(diffs);
        if min_diff <= tolerance && ~matched(idx)
            correct = correct + 1;
            matched(idx) = true;  % Prevent double matching
        end
    end

    accuracy = (correct / length(ref_peaks)) * 100;
end

