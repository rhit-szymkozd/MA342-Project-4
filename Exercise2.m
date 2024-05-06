clc
clear variables
close all

% stocks = hist_stock_data('30082023','30042024','AAPL');
% AAPL = process_stock_data(stocks,0);
% save("AAPL.mat","AAPL")

load AAPL.mat
AAPL = AAPL(end-19:end);
% APPL = AAPL.Open(end-20:end,:);
[mu, sigma] = genMuSigma(AAPL);

%% Part 1
% P1 = generateStock12_2(AAPL, 20, 1);
load P1_AAPL_12_2.mat
[mu1, sigma1] = genMuSigma(P1);

%% Part 2
% P2 = generateStock12_3(AAPL, 20, 1);
load P1_AAPL_12_3.mat
[mu2, sigma2] = genMuSigma(P2);

fprintf("original data: mu = %0.4f, sigma = %0.4f\n", mu, sigma);
fprintf("P1: mu = %0.4f, sigma = %0.4f\n", mu1, sigma1);
fprintf("P2: mu = %0.4f, sigma = %0.4f\n", mu2, sigma2);

%% Part 3
% variance estimated with 12.3 is over 3x larger than it should be. 
% 12.3 has variance increase with time, exponentially affects estimate
% remove sqrt(t) term so that variance will not increase with time

% P3 = generateStock12_3_2(AAPL, 20, 1);
load P1_AAPL_12_3v2.mat
[mu3, sigma3] = genMuSigma(P3);

fprintf("P3: mu = %0.4f, sigma = %0.4f\n", mu3, sigma3);


%% FUNCTIONS

function [mu, sigma] = genMuSigma(P)

R = (P(2:end,:)-P(1:end-1,:))./P(1:end-1,:);
mu = mean(R);
sigma = std(R);

end

function [P] = generateStock12_2(P, T, plott)
histLength = length(P);
[mu, sigma] = genMuSigma(P);
deltat = 1;
for i = histLength:histLength+T-1
    phi = normrnd(0, 1);
    P(i + 1) = P(i) + mu*P(i)*deltat + sigma*P(i)*sqrt(deltat)*phi;
end
if(plott)
    figure()
    hold on
    plot(1:histLength,P(1:histLength),'b-');
    plot(histLength:length(P),P(histLength:end),'r-');
    xlabel("Time (Days)");
    ylabel("Stock Price ($USD")
    legend("Historical Price (Open)", "Simulated Price","Location", "northwest")
end
end

function [P] = generateStock12_3(P, T, plott)
histLength = length(P);
[mu, sigma] = genMuSigma(P);
t=1:T;
phi = normrnd(0, 1, 1, T);
append = P(end).*exp((mu-sigma^2/2).*t + sigma.*sqrt(t).*phi);
P = [P;append'];

if(plott)
    figure()
    hold on
    plot(1:histLength,P(1:histLength),'b-');
    plot(histLength:length(P),P(histLength:end),'r-');
    xlabel("Time (Days)");
    ylabel("Stock Price ($USD")
    legend("Historical Price (Open)", "Simulated Price","Location", "northwest")
end

end

function [P] = generateStock12_3_2(P, T, plott)
histLength = length(P);
[mu, sigma] = genMuSigma(P);
t=1:1:T;
phi = normrnd(0, 1, 1, length(t));
append = P(end).*exp((mu-sigma^2/2).*t + sigma.*phi);
P = [P;append'];

if(plott)
    figure()
    hold on
    plot(1:histLength,P(1:histLength),'b-');
    plot([histLength, t+histLength],P(histLength:end),'r-');
    xlabel("Time (Days)");
    ylabel("Stock Price ($USD")
    legend("Historical Price (Open)", "Simulated Price","Location", "northwest")
end

end