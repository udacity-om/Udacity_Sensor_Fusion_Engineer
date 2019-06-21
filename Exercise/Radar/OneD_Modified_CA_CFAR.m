% Implement 1D CFAR using lagging cells on the given noise and target scenario.

% Close and delete all currently open figures
close all;
clearvars;
clc;

% Data_points
Ns = 1000;

% Generate random noise
s=abs(randn(Ns,1));

%Targets location. Assigning bin 100, 200, 300 and 700 as Targets with the amplitudes of 8, 9, 4, 11.
s([100 ,200, 300, 700])=[8 9 4 11];

%plot the output
plot(s);

% TODO: Apply CFAR to detect the targets by filtering the noise.

% 1. Define the following:
% 1a. Training Cells
T = 12;
% 1b. Guard Cells 
G = 4;

% Offset : Adding room above noise threshold for desired SNR 
offset=5;

% Vector to hold threshold values 
threshold_cfar = [];

%Vector to hold final signal after thresholding
signal_cfar = [];

% 2. Slide window across the signal length
for i = 1:(Ns-(G+T))     

    % 2. - 5. Determine the noise threshold by measuring it within the training cells
    if(i <= (G+1))
        lagging_average = 0;
    else
        lagging_training_cell_min_index = max(1, i-G-T);
        lagging_average = mean(s(lagging_training_cell_min_index:i-G-1));
    end
    
    leading_average = mean(s(i+G+1:i+G+T));
    
    %average_noise = mean([lagging_average leading_average]);
    average_noise = max(lagging_average, leading_average);
    noise_threshold = average_noise * offset;
    threshold_cfar = [threshold_cfar, {noise_threshold}];

    % 6. Measuring the signal within the CUT
    CUT = s(i);

    % 8. Filter the signal above the threshold
    if(CUT > noise_threshold)
        signal = CUT;
    else
        signal = 0;
    end  
    signal_cfar = [signal_cfar, {signal}];
end


% plot the filtered signal
plot (cell2mat(signal_cfar),'g--');

% plot original sig, threshold and filtered signal within the same figure.
figure,plot(s);
hold on,plot(cell2mat(circshift(threshold_cfar,G)),'r--','LineWidth',2)
hold on, plot (cell2mat(circshift(signal_cfar,(T+G))),'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')