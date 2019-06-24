clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
% Velocity Resolution = 0.5 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User Defined Range and Velocity of target

R=110;     %initial distance of the target
v=-30;      %speed of the target 

%% FMCW Waveform Generation

%Designing the FMCW waveform by giving the specs of each of its parameters.

%Maximum Range supported by the Radar. This will in turn define each chirp
%sweep time. Here we limit the range of radar to 200m considering the
%limited power on RF transmitter and the path loss. 
R_max = 200;

%range_resolution
delta_r = 1;  %m

%Speed of light
c=3e8;                   

%Sweep time for each chirp is defined as rule by 5.5 times of round trip
%time. Sweep time = 5.5*2*Rmax/c
Swp_Tm = 5.5*2*R_max/c; 

%Sampling frequency mainly used here for plotting the chirps
Fs=150e6;                

%Operating carrier frequency of Radar - 77GHz. The ADAS frequency of operation for
%radar is from 76GHz to 81GHz. 
fc=77e9;                 %carrier freq

%Bandwidth for the each chirp is defined by the range resolution spec
%the radar. range_resolution = c/2B. Range resolution is the capability of
%the radar to resolve two targets in range. To resolve two targets 1 m
%apart the BW = c/2*range_resolution = 150 MHz
BW=c/2*delta_r;                %sweep freq (bandwidth) for each chirp

%The slope of the chirp given by BW/Swp_Tm
Slope=BW/Swp_Tm;         %sweep rate (slope)
disp(Slope);

%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. Also, called Doppler bins or Slower samples
D=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. This gives the Range FFT bins. Also,
%called faster samples
N=1024;                  %for length of time OR # of range cells

% The total time to send all the D chirps with N sampling on each
t=linspace(0,D*Swp_Tm,D*N); %total time for samples



%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t));

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Spectrum Analyzer
%Enabling the Spectrum Analyzer (function provided directly to students)

specanalyzer = dsp.SpectrumAnalyzer('SampleRate',Fs,...
    'PlotAsTwoSidedSpectrum',false,...
    'Title','Spectrum for received and dechirped signal',...
    'ShowLegend',true);

%% Radar Trasmission and Reception
% Running the radar scenario over the time. 

for i=1:length(t)
    
    %For each time sample update the Range and Velocity of the Target. 
    % Also, update the return time for RF signal to reach to the target and
    % be received by the Rx. Rnew = R + v(t) and return time update would
    % be td = 2*Rnew/c
    
    r_t(i)=R+(v*t(i)); % range of the target in terms of its velocity and initial range
    td(i)= 2*r_t(i)/c; % delay for received signal
    
    %Besides, for each time sample we need to update the transmitted and
    %received signal. Since these are time domain signal they will vary in
    %amplitude and phase with time. 
    
    % Transmitted Signal  = SIN(2pi(fc*t + Slope*(t^2)/2)
    % Received Signal - using the delayed version of the transmit waveform
    % xt = x(t-tau);
    
    Tx(i)=cos(2*pi*(fc*t(i)+.5*Slope*t(i)^2)); %transmitted signal
    Rx(i)=cos(2*pi*(fc*(t(i)-td(i))+.5*Slope*(t(i)-td(i))^2)); %received signal
    
    
    %Dechirping of the signal - this is where we calculate the shift between
    %the transmit and beat frequency. The shift gives us the range and
    %doppler information. Also, called the demixing of the signal.
    %To get the dechirped signal we need element by element multiplication
    %of Tx and Rx waveforms.
    
    %This would give both doppler and range post FFT
    
    Mix(i)=(Tx(i).*Rx(i));
    
    %Plotting the output of the dechirped signal (IF signal).
%     specanalyzer(Mix(i));

end

%% plotting the chirps sequences from Tx, Rx and Mix signal. 

figure ('Name','Transmitted, Received and Beat Signals');
subplot(3,2,1);
spectrogram(Tx,256,250,256,Fs,'yaxis');
subplot(3,2,2);
spectrogram(Mix,256,250,256,Fs,'yaxis')
subplot(3,2,4);
spectrogram(Rx,256,250,256,Fs,'yaxis');

%% RANGE MEASUREMENT

%reshape the vector into N*D array. N and D here would define the size of
%Range and Doppler FFT respectively.

Mix=reshape(Mix,[length(Mix)/D,D]);

%Determing the FFT size and it has to be closest exponent of 2.
[My,Ny]=size(Mix');
nfft = 2^nextpow2(N);

% To determing Range perform the 1st FFT across non-doppler dimension. 
% The FFT is performed to convert the time domain signal into frequency
% domain.

Mix1 = fft(Mix,nfft)/Ny;

% The output of FFT is in complex domain so it has magnitude and phase information.
% Here we use the magnitude so we take the absolute value to plot the range
Mix1 = abs(Mix1);


% Output of FFT is double sided singal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
Mix1 = Mix1(1:nfft/2);

%Since we are using just one half we multiply the values by 2.
Mix1(2:end-1) = 2*Mix1(2:end-1);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)
plot(Mix1);
axis ([0 500 0 1]);



%% RANGE DOPPLER RESPONSE


% 2 Dimensional FFT is needed to get the Range Doppler response. As seen
% above the 1st dimension gives the range along 1st dimension of the matrix
% and the 2nd FFT on the other dimension gives the velocity.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.
doppler_axis = linspace(-.5,0.5-1/My,My)*2*100;
range_axis = linspace(-200,200,nfft/2)*((nfft/2)/400);

% 2D FFT using the FFT size for both dimensions.
Y = fft2(Mix,nfft,My);

% Taking just one side of signal from Range dimension.
Y = Y(1:nfft/2,1:My);
Y = fftshift (Y);
Y = abs(Y);
Y = 10*log10(Y);


%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
figure,surf(doppler_axis,range_axis,Y);


%% CFAR implementation

%In the radar receiver, the returning echoes are typically received by the 
%antenna, amplified, down-converted and then passed through detector circuitry 
%that extracts the envelope of the signal (known as the video signal). 
%This video signal is proportional to the power of the received echo and comprises
%the wanted echo signal and the unwanted power from internal receiver noise 
%and external clutter and interference.

%The role of the constant false alarm rate circuitry is to determine the 
%power threshold above which any return can be considered to probably originate from a target. 
%If this threshold is too low, then more targets will be detected at 
%the expense of increased numbers of false alarms. 
%Conversely, if the threshold is too high, then fewer targets will be detected, 
%but the number of false alarms will also be low. 
%In most radar detectors, the threshold is set in order to achieve a required probability of false alarm 
%(or equivalently, false alarm rate or time between false alarms).

% disp (max(Y));

%Slide Window through the complete Range Doppler Map
Nr = nfft/2;
Nd = My;

%Select the number of Training Cells in Range as well as Doppler dimension 
%to run the Noise estimate
Tr=16;
Td=12;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr=4;
Gd=4;

%Offset the threshold by 
offset=6;

threshold_cfar = [];

%To store noise_level for each iteration over training cell
noise_level = zeros(1,1);

%Using the Range doppler output from the 2DFFT
signal = Y; 

%Set high value for above threshold cells
max_T = 60;

%slide window across all the bins of Range Doppler Map (RDM)
for i = 1:(Nr-(Gr+Tr+1))
    for j = 1:(Nd-(Gd+Td+1))
        %Measure the noise across the Training Cells
        for m = i:i+Tr-1
            for q = j:j+Td-1
                noise_level=noise_level + db2pow(Y(m,q)); %for addition convert the dB values into pow
%                 disp(noise_level);
            end
        end
        %Calculate the threshold based on the noise level measured through
        %training cells and adding an optimal offset.
        threshold = pow2db(noise_level/(Td*Tr))+offset; % After addition revert it back into dB
%         threshold_cfar = [threshold_cfar, {threshold}];
        
        %Measure the signal in Cell Under Test (CUT)and compare against
        %Threshold. If the Signal is above the threshold then set its value
        %to high, if not then set it to 0.
        signal_T=Y(i+Tr+Gr+1,j+Td+Gd+1);
        if(signal_T<threshold)
            signal(i+Tr+Gr+1,j+Td+Gd+1)=0;
        else 
            signal(i+Td+Gd+1,j+Tr+Gr+1)=max_T;
        end
        noise_level = 0;
    end
end

% The CUT cell matrix will be smaller than the Range Doppler Map, so it
% few cells will not be thresholded. To keep the map size same set those
% values to 0. 

for a = 1:Nr
    for b = 1:Nd
        if (signal(a,b)~=0 && signal(a,b)~=max_T)
            signal(a,b) = 0;
        end
    end
end
% 
%  disp(max(signal));

%disp(numel(find((signal==max_T))));

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,signal);
colorbar;

%% Clustering

%Here merge multiple detections suspected to be of the same vehicle to a single detection. Look for
%detections that are closer than the size of a vehicle. 
%Detections that fit this criterion are considered a cluster and are merged
%to a single detection 

%Determine all the detections and their range. Based on the difference of
%their range club them together under the same group.


%Find the detections using the 'find' function.
 detections = find(signal==max_T);
 
 %Determine the row and column of each detection.
 [row,col] = find(signal==max_T);
 
 %# of detections
 num_detections = length(row);
 
 %determine the range based on range/doppler bin estimate for each range
 %and doppler values
 
 if (row>256)
     range = (row-256)*1;
 else 
     range = -(row)*1;
 end
 if (col>64)
     velocity = (col-64)*100/64;
 else
     velocity = -(64-col)*100/64;
 end
 
 
 %generate the distance_delta vector
 distance_delta = zeros(length(row));

 %length of vehicle
 vehicle_length = 2.5; %m
 
 % measure the distance between all sets of the detections.
  for m = 1:num_detections
     cluster_ids(m) = m;
     for n = 1+m:num_detections
         distance_delta(m,n) = range(m)- range(n);
         
         % if distance is less than the length of vehicle group them
         % together 
         if (distance_delta(m,n)< vehicle_length) 
               cluster_ids(m)=1;
         else
               cluster_ids(m)=m;
         end
    end
  end
 
  %disp(cluster_ids);