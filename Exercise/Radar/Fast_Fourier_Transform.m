clearvars;
clc;

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

% TODO: Form a signal containing a 77 Hz sinusoid of amplitude 0.7 and a 43Hz sinusoid of amplitude 2.
f1 = 50;
f2 = 120;
A1 = 0.7;
A2 = 2;
S = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t);

% Plot the pure signal in the time domain.
subplot(221);
plot(1000*t(1:50) ,S(1:50))
title('Signal')
xlabel('t (milliseconds)')
ylabel('S(t)')

% Run fft
S_fft = fft_func(S, L);

% Plotting. The amplitudes are exactly at 0.7 and 2.
f = Fs*(0:(L/2))/L;
subplot(222);
plot(f,S_fft) 
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|Sfft(f)|')

% Corrupt the signal with noise 
X = S + 2*randn(size(t));

% Plot the noisy signal in the time domain. It is difficult to identify 
% the frequency components by looking at the signal X(t). 
subplot(223);
plot(1000*t(1:50) ,X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

% Run fft
X_fft = fft_func(X, L);

% Plotting. The amplitudes are not exactly at 0.7 and 1, as expected, because of 
% the added noise. On average, longer signals produce better frequency 
% approximations.
f = Fs*(0:(L/2))/L;
subplot(224);
plot(f,X_fft) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|Xfft(f)|')

function X = fft_func(x, L)
    % Run the fft for the signal using MATLAB fft function for dimension of samples N
    X = fft(x);

    % The output of FFT processing of a signal is a complex number (a+jb). 
    % Since, we just care about the magnitude we take the absolute value (sqrt(a^2+b^2)) of the complex number.
    X = abs(X/L);

    % FFT output generates a mirror image of the signal. But we are only 
    % interested in the positive half of signal length L, since it is the 
    % replica of negative half and has all the information we need.
    X = X(1:L/2+1);
    X(2:end-1) = 2*X(2:end-1);
end

