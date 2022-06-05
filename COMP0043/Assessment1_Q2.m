%% COMP0043 - Numerical Methods for Finance

%  Date: 13th December 2021
%  Assessment 1: Online Remote Test (40%)

% The following is the script for QUESTION 2. A separate file has been
% provided for Question 1.

clear all
close all

% QUESTION 2

% Option parameters
T = 0.5 ; % Time period
K = 1 ; % Strike price
r = 0.1 ; % Risk free interest rate
q = 0 ; % Dividend rate
S0 = 1 ; % Initial Stock Price
N = 50 ; % Monitoring dates

% GBM Process Parameters
sigma = 0.3 ; % Volatility/Difffusion
dt = T/N ; % Size of our timesteps
t = 0:dt:T ; % Discretization of our grid

% Risk Neutral Measure
muRN = r-q-0.5*sigma^2 ;

% Monte Carlo Parameters
nblocks = 5000 ;
nsample = 5000 ;

%% Question 2 -- Part (a) 
% Monte Carlo Simulation - Lookback Call Option

for i = [1:nblocks]
    % We compute our traditional ABM on log-price
    dX = muRN*dt + sigma*sqrt(dt)*randn([nsample,N]) ;
    
    % Now we need to cumulatively sum the values over the time steps to get
    % each path
    X = [zeros([nsample,1]) , cumsum(dX,2)] ;  

    % We now take its exponential to recover the stock price
    S = S0*exp(X) ;

    % Compute the payoff of the option at each of the timesteps (monitoring
    % points)
    S2 = max(S-K,0) ;
    
    % Take the maximum payoff over the entire monitoring period
    for j = [1:nsample]
        S3(j) = max(S2(j,:)) ;
    end

    % Compute the discounted payoff for that path
    VCMCb(i) = exp(-r*T)*mean(S3) ;
end

VCMC = mean(VCMCb) ; 

%% Question 2 -- Part (b)
% Monte Carlo Simulation - Lookback Put Option

for i = [1:nblocks]
    % We compute our traditional ABM on log-price
    dX = muRN*dt + sigma*sqrt(dt)*randn([nsample,N]) ;
    
    % Now we need to cumulatively sum the values over the time steps to get
    % each path
    X = [zeros([nsample,1]) , cumsum(dX,2)] ;  

    % We now take its exponential to recover the stock price
    S = S0*exp(X) ;

    % Compute the payoff of the option at each of the timesteps (monitoring
    % points)
    S2 = max(K-S,0) ;
    
    % Take the maximum payoff over the entire monitoring period
    for j = [1:nsample]
        S3(j) = max(S2(j,:)) ;
    end

    % Compute the discounted payoff for that path
    VPMCb(i) = exp(-r*T)*mean(S3) ;
end

VPMC = mean(VPMCb) ; 



