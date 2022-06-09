%% COMP0043 - Numerical Methods for Finance

%  Date: 13th December 2021
%  Assessment 1: Online Remote Test (40%)

% The following is the script for QUESTION 1. A separate file has been
% provided for Question 2.

clear all
close all

% QUESTION 1 -- Part (a)

% Parameters
d = 5 ;
lambda = 5 ;
nsample = 10^5 ;
deltax = 0.2 ;

% Define the size of the grid (and therefore the bins)
x = 0:deltax:20 ;

% PDF for the N.C. Chi-Sq. Distribution using MATLAB function
f = pdf('ncx2',x,d,lambda) ;

% Sampling the N.C. Chi-Sq. using MATLAB function
R = ncx2rnd(d,lambda,[1,nsample]) ;

% Plot the normalized histogram
figure(1)
histogram(R,'BinEdges',x,'Normalization','pdf')
hold on
plot(x,f,'r')
xlim([0,20])
xlabel('x')
ylabel('f')
legend('Sampled','Theory')
title('N.C. Chi-Sq. Distribution with d=5 and \lambda=5')


% QUESTION 1 -- Part (b)

% Generate standard uniform random numbers
U1 = rand([1000,1]) ;

% Plug these into the iCDF function to obtain N.C. Chi-Sq. random numbers
% with degrees of freedom and lambda as required
X = icdf('ncx2',U1,d,lambda) ;

% Take those N.C. Chi-Sq. random numbers and put them into the PDF 
fX = ncx2pdf(X,d,lambda) ;

U2 = rand([1000,1]) ;

% Multiply the PDF numbers generated by standard unifrom random numbers to
% scatter them evenly under the PDF, as required
X2 = U2.*fX ; 

figure(2)
plot(X,X2,'b.')
xlim([0,20])
xlabel('x')
ylabel('f')
title('N.C. Chi-Sq. Distributed Random Numbers with d=5 and \lambda=5')






