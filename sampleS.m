function [B, S] = sampleS(a, b, dt, S0, n)
%
% function S = sampleS(a, b, T, dt, S0, n)
% 
% A function that simulates stock prices that follow geometric Brownian
% motion, given parameters.
%
% The function uses Euler's simulation scheme
% 
% a - drift parameter
% b - volatility parameter
% dt - mesh size
% S0 - intial stock price
% n - number of simulated points




B = zeros(n,1);     % initialize Brownian motion
S = zeros(n,1);     % initialize Asion option price
S(1) = S0;          % initial price

B(2:end) = cumsum(sqrt(dt)*randn(n-1,1));       % Brownian motion path
S = S0*exp(a*dt*(0:n-1)' + b*B);                % stock price spath