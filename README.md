# ECG-Signal-Denoising-and-Analysis-using-Static-and-Adaptive-Filtering
This project focuses on ECG signal denoising and analysis using the MIT-BIH Arrhythmia Database. The goal is to simulate common ECG noises, apply static filters (IIR elliptic) and adaptive filters (LMS), and compare their performance using SNR, RMSE, and R-peak detection accuracy.

# Key Features
Noise Simulation: Adds realistic ECG noise types:
Powerline interference (50 Hz)
Baseline wander (0.3 Hz)
EMG/muscle artifact noise
Frequency Analysis: Uses FFT, spectrograms, and STFT for spectral analysis.
Static Filtering: Applies notch, high-pass, and low-pass elliptic filters.
Adaptive Filtering: Implements LMS adaptive filtering with reference signals for powerline and baseline wander removal.
Quantitative Metrics: Computes SNR, RMSE, and R-peak detection accuracy for clean, noisy, and filtered signals.
Visualization: Plots signals at every stage (clean, noisy, filtered) along with frequency and time-frequency analysis.

# Requirements
MATLAB (R2020 or later recommended)
Signal Processing Toolbox
DSP System Toolbox
MIT-BIH Arrhythmia Database (WFDB Toolbox installed)

# Usage
Set your dataset directory in the script:
cd('D:\path\to\mit-bih-arrhythmia-database-1.0.0');  
Run the script in MATLAB.

# The program will:
Load and slice a 3-second ECG segment
Add synthetic noises
Apply filtering (static & adaptive)
Plot time/frequency domain results
Display SNR, RMSE, and R-peak detection accuracy

# Outputs
Time-domain plots: clean, noisy, filtered signals
Frequency spectra and spectrograms
STFT-based time-frequency plots
Comparison of static vs. adaptive filtering

# Performance metrics:
Initial vs. static vs. adaptive SNR
RMSE values
R-peak detection accuracy (%)

# Example Results
Static Filters: Effective at reducing powerline noise but less adaptive to baseline drift.
Adaptive Filters: More robust against multiple noise sources and achieve higher accuracy in R-peak detection.

# Future Improvements
Real-time filtering implementation
Application to longer ECG signals (full 30 min)

Comparison with wavelet-based denoising

Integration with QRS detection algorithms
