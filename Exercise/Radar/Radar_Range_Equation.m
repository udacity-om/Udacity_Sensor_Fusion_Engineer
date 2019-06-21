% Operating frequency (Hz)
fc = 77.0e9;

% Transmitted power (W)
Pt = 3e-3;

% Antenna Gain (linear)
G = 10000;

% Minimum Detectable Power
Ps = 1e-10;

% RCS of a car
RCS = 100;

%Speed of light
c = 3*10^8;

% TODO: Calulate the wavelength
lambda = c/fc;

% TODO: Measure the Maximum Range a radar can see
R = power(((Pt*power(G, 2)*power(lambda, 2)*RCS)/(Ps*power((4*pi), 3))), 1/4);

disp(R)