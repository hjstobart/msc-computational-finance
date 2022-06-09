%% COMP0043 - Numerical Methods for Finance

%  Date: 11th May 2022
%  Assessment 2: Summer Exam (60%)

% ===========================
%     --- Question 3 ---
% ===========================

%% Define the grids

clear all
close all

% GRID IN REAL SPACE
N = 2048; % The number of grid points
dx = 0.01; % Step size of the grid in real space
upperx = N*dx; % Upper truncation limit in real space
x = dx*(-N/2:N/2-1); % Grid in real space

% GRID IN FOURIER SPACE (Pulsation)
dxi = (2*pi)/(N*dx); % Step size of the grid in fourier space
upperxi = N*dxi; % Upper truncation limit in fourier space
xi = dxi*(-N/2:N/2-1); % Grid in fourier space

%% FFT Implementation (from script exp2lor.m)

tic

% ANALYTICAL expressions
% -----------------------------
a = 1; % Activity parameter

fa = 0.5*a*exp(-a*abs(x)); % Laplace
% We will use this to check that the inverse numerical FFT does a good
% approximation of the analytical expression

Fa = a^2./(a^2 + xi.^2); % Lorentz (Pulsation)

% NUMERICAL approximations
% -----------------------------
Fn = fftshift(ifft(ifftshift(fa)))*upperx;
fn = fftshift(fft(ifftshift(Fa)))/upperx;

cputime_fft = toc;

%% Discrete Fourier Transform

% The discrete Fourier Transform takes the analytical expression of
%    f(x): Bilateral Exponential / Laplacian
% ...and computes...
%    f(xi): Lorentzian / Cauchy

tic

% Initialize vector
F_DFT = zeros([1,N]);

for k = 1:N
    F_DFT(k) = sum( exp(1i*xi(k)*x) .* (fa * dx));
end

cputime_DFT = toc;

%% Printing the results to screen

fprintf('%20s%15s%15s\n','','FFT','DFT');
for j = 1:5
    fprintf('%20s%15.12f%15.12f\n','',real(Fn(j)),real(F_DFT(j)))
end
fprintf('%20s%15.8f%15.8f\n\n','CPU Time (s)',cputime_fft,cputime_DFT)


%% Figure

figure(1), clf, 
hold on;
plot(xi,real(Fn),'b')
plot(xi,imag(Fn),'m')
plot(xi,Fa,'c:')
plot(xi,real(F_DFT),'ko')
plot(xi,imag(F_DFT),'go')
axis([-10 10 0 1.1])
title('Analytical FT, FFT and DFT of the Double Exponential function',Interpreter='latex')
xlabel('\xi')
ylabel('f(\xi)')
legend('Re(Fn)','Im(Fn)','Fa','Re(DFT)','Im(DFT)')

















