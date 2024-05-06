clc
clear variables
close all

% stocks = hist_stock_data('01072014','01122014','groupTickers.txt');
% save("chosenStocks.mat","stocks")
load("chosenStocks.mat");
P = process_stock_data(stocks,0);

[numDays, numStocks] = size(P);
R = (P(2:end,:)-P(1:end-1,:))./P(1:end-1,:);
r = mean(R);
C = cov(R);

a = (0:0.01:1).^2;
for i = 1:length(a)
alpha = a(i);
[w, optVal] = quadprog((1-alpha)*2*C, -alpha*r, [], [], ones(1,numStocks), 1, zeros(numStocks,1), ones(numStocks,1));
[w2, optVal2] = quadprog((1-alpha)*2*C, -alpha*r, [], [], ones(1,numStocks), 1, zeros(numStocks,1), 0.2*ones(numStocks,1));
risk(i) = w'*C*w;
ret(i) = r*w;

risk2(i) = w2'*C*w2;
ret2(i) = r*w2;
end
hold on
plot(risk, ret, 'b-');
plot(risk2, ret2, 'r-');
