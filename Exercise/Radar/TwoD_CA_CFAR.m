% Implement 1D CFAR using lagging cells on the given noise and target scenario.

% Close and delete all currently open figures
close all;
clearvars;
clc;

% Data_points
R = 100; % In range dimension
D = 100; % in doppler dimension

% Generate random noise
s=abs(randn(R,D));

%Targets location. Assigning bin 100, 200, 300 and 700 as Targets with the amplitudes of 8, 9, 4, 11.
s(10, 10) = 8;
s(20, 20) = 9;
s(30, 30) = 4;
s(70, 70) = 11;

%plot the output
imagesc(s);

% TODO: Apply CFAR to detect the targets by filtering the noise.

% 1. Define the following:
% 1a. Training Cells
Tr = 3;
Td = 3;
% 1b. Guard Cells 
Gr = 1;
Gd = 1;

% Offset : Adding room above noise threshold for desired SNR 
offset=6;

% 2. Slide window across the signal length
cfarWin = ones(((Td+Gd)*2 + 1), (Tr+Gr)*2 + 1);
cfarWin(Td+1:end-Td, Tr+1:end-Tr) = 0;
cfarWin = cfarWin/sum(cfarWin, 'all');

noise_level = conv2(s, cfarWin, 'same');
threshold_cfar = noise_level*offset;

detection = (s > threshold_cfar);
detection = detection.*s;

% plot the filtered signal
figure, imagesc(detection);

% plot original sig, threshold and filtered signal within the same figure.
figure;
surf(s)
hold on
surf(threshold_cfar)
surf(detection)