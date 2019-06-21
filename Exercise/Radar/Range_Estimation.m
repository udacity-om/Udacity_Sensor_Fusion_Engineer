% TODO : Find the Bsweep of chirp for 1 m resolution
d_res = 1;
B_sweep = c/(2*d_res);

% TODO : Calculate the chirp time based on the Radar's Max Range
R_max = 300;
c = 3.0e8;
Ts = (2*R_max)/c;

% TODO : define the frequency shifts 
fb = [0, 1.1, 13, 24]; % MHz
calculated_range = (c*Ts*fb)/(2*B_sweep);

% Display the calculated range
disp(calculated_range);