% 2-D Transform
% The 2-D Fourier transform is useful for processing 2-D signals and other 2-D data such as images.
% Create and plot 2-D data with repeated blocks.

P = peaks(20);
X = repmat(P,[5 10]);
subplot(121);
imagesc(X)

% TODO : Compute the 2-D Fourier transform of the data.  
X_fft = fft2(X);

% Shift the zero-frequency component to the center of the output, and 
X_fft = fftshift(X_fft);
X_fft = abs(X_fft);

% plot the resulting 100-by-200 matrix, which is the same size as X.
subplot(122);
imagesc(X_fft);