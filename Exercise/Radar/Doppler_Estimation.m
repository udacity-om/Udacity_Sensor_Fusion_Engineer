% Doppler Velocity Calculation
c = 3*10^8;         %speed of light
frequency = 77e9;   %frequency in Hz

% TODO : Calculate the wavelength
lambda = c/frequency;

% TODO : Define the doppler shifts in Hz using the information from above 
fd = [3, 4.5, 11, -3]; %KHz

% TODO : Calculate the velocity of the targets  fd = 2*vr/lambda
vr = (fd*lambda)/2;

% TODO: Display results
disp(vr);