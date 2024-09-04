function [A,B] = fourierCoeff(x)
% Returns the array A and B which contain the trigonometric fourier
% coefficents of the fourier series representation of signal x

L = length(x);
fft_x = fft(x,L-1);

% Normalized fft
nfft_x = fft_x/L;

% cutoff defined to discard redundant half of FFT
cutoff = ceil(L/2);

A = 2*real(nfft_x(1:cutoff));
B = -2*imag(nfft_x(1:cutoff));

% a0 in the fourier series
A(1) = A(1)/2;

end