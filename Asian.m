clear, clc, close all;

% Asian options parameters - parameters of the Brownian motion
a = -.045; b = .3;
S0 = 8;

% Simulation parameters
T = 100;
dt = (1/.3*log(.01))^2 / 1000;
% dt = .1^2/100;
n = ceil(T/dt);

fig = figure;
col = ['bgrmy'];


for i=1:5
    % Sample Brownian and geometric Brownian motion
    [B, S] = sampleS(-.045, .3, dt, S0, n);
    
    % Plot
    x_plot = round(linspace(1, n, 1000));   % x axis
    
    subplot(211);
    title('Standard Brownian Motion');
    plot(dt*x_plot, B(x_plot), col(i)); hold on;
    xlabel('Time'); ylabel('B(t)');
    
    subplot(212);
    title('Geometric Brownian Motion');
    plot(dt*x_plot, S(x_plot), col(i)); hold on;
    xlabel('Time'); ylabel('S(t)');
end



%% Asian option with strike price 10, maturity date T=30;

K = 10;         % strike price
T = 30;         % maturity date
S0 = 8;         % initial stock price

dt = (1/.3*log(1.05))^2;        % simulation step size (mesh size)
n = ceil(T/dt);                 % number of simulated points

R = 5e6;                        % number of sample path simulations
X = zeros(R, 1);                % initialize option price (at time T)


% drawing Brownian motion samples (sample paths) and integrating asian
% option price using the trapezoidal rule
for i=1:R
    [B, S] = sampleS(-.045, .3, dt, S0, n);
    X(i) = max(0, 1/T*trapz(dt*(0:n-1), S) - K);
end


fprintf('Asian option price: %3.4f    ', mean(X))
fprintf('Confidence interval <%3.4f, %3.4f>\n', ...
    (mean(X) - 1.96*std(X)/sqrt(R)), (mean(X) + 1.96*std(X)/sqrt(R)))
