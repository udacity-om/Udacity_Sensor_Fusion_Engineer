close all;
clearvars;
clc;

max_T = 10e3;
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
c = 3e8;
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
 R = 110;
 v = -20;

%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
R_res = 1;
R_max = 200;

B = c/(2*R_res);
Tchirp = 5.5*(2*R_max)/c;
slope = B/Tchirp; % For given system requirements the calculated slope should be around 2e13

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)  
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    current_time = t(i);
    if(i == 1)
        elapsed_time = current_time;
    else
        elapsed_time = current_time - t(i-1);
    end
    
    R = R + v*elapsed_time;
    
    % trip time is the time it took the Tx signal to reach the target(initailly at
    % 110m) and get reflected back to sensor
    trip_time = 2*R/c;
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*((fc*current_time) + ((slope*power(current_time,2))/2)));
    time_diff = current_time-trip_time;
    Rx(i) = cos(2*pi*((fc*time_diff) + ((slope*power(time_diff,2))/2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i)*Rx(i);
    
end

%% RANGE MEASUREMENT
 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr,Nd]);
%plotting the received signal
figure ('Name','Received signal in time domain')
plot(Mix(:,1));
title('Received signal in time domain');
xlabel('Time');
ylabel('Amplitude');


 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
first_fft = fft(Mix, Nr, 1);

 % *%TODO* :
% Take the absolute value of FFT output
first_fft = abs(first_fft/Nr);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
first_fft = first_fft(1:Nr/2+1, :);
first_fft(2:end-1, :) = 2*first_fft(2:end-1, :);

%plotting the range
figure ('Name','Range from First FFT')

% *%TODO* :
% plot FFT output 
plot(first_fft(:,2));
axis ([0 200 0 1]);
title('Range from First FFT');
xlabel('Range');
ylabel('Amplitude');

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure('Name','Range-Doppler from 2D FFT');
surf(doppler_axis,range_axis,RDM);
title('Range-Doppler from 2D FFT');
xlabel('Doppler');
ylabel('Range');
zlabel('Amplitude');

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 3;
Gd = 6; % As the signal leaks more in doppler direction

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 2; % After the guard cells the noise is almost uniform. So a small region is enough to find the mean noise
Td = 2;

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 15; % Looked at the difference in peak and noise value and it comes up to ~17dB

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR

for i = Tr+Gr+1:(Nr/2)-(Gr+Tr)
   for j = Td+Gd+1:Nd-(Gd+Td)
       noise_level = 0;

       for p = i-(Tr+Gr):i+(Tr+Gr)
           for q = j-(Td+Gd):j+(Td+Gd)
               if((abs(i-p) > Gr) || (abs(j-q) > Gd))
                   noise_level = noise_level + db2pow(RDM(p,q));
               end
           end
       end

       threshold = pow2db(noise_level/((2*(Td+Gd)+1)*(2*(Tr+Gr)+1) - (2*Gr + 1)*(2*Gd + 1)));
       threshold = threshold + offset;

       CUT = RDM(i,j);
       if(CUT < threshold)
           RDM(i,j) = 0;
       else
           RDM(i,j) = max_T;
       end           
   end
end

% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 RDM(RDM ~= 0 & RDM ~= max_T) = 0;
 RDM(RDM == max_T) = 1;

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure('Name','Detection after 2D-CFAR');
surf(doppler_axis,range_axis,RDM);
title('Detection after 2D-CFAR');
xlabel('Doppler');
ylabel('Range');
zlabel('Amplitude');
colorbar;


 
 