clc
clear all
close all

start_date = '25102023';
end_date = '25042024';
varargin = {"MSFT","AAPL","AMZN","JPM","V","WMT","UNH","PG","JNJ","HD","MRK","CVX","CRM","KO","DIS","MCD","CSCO","CAT","AXP","IBM","VZ","AMGN","INTC","NKE","GS","HON","BA","MMM","TRV","DOW"};
data_raw = hist_stock_data(start_date,end_date,varargin);

n = length(cell2mat(data_raw(3,1,1)));

data = zeros(n,length(varargin));

for(i=1:length(varargin))
    data(:,i) = process_stock_data(data_raw(:,i),0);
end

R = (data(2:end,:)-data(1:end-1,:))./data(1:end-1,:);%returns

r = mean(R)';

C = cov(R);

n = length(C);

alpha = [0:0.01:1].^2;

mult = [1,0.5,0.2];

for(multi = 1:3)
for(j=1:length(alpha))

    [w(:,j),optVal] = quadprog((1-alpha(j)).*2.*C,-alpha(j).*r,[],[],ones(1,n),[1],zeros(n,1),mult(multi).*ones(n,1));
    risk(j) = w(:,j)'*C*w(:,j);
    expected_return(j) = r'*w(:,j);

end

slope = (expected_return(2:end)-expected_return(1:end-1))./(risk(2:end)-risk(1:end-1));

[lowest,lowest_i] = min(abs(slope-1));%find where slope is closest to 1
lowest = lowest+1
best_alpha = alpha(lowest_i)%find slope and alpha values at the best point
best_w = w(:,lowest_i)%best optimal portfolio
count = 1;
for(k=1:length(best_w))
if(best_w(k)>0.0001)%find any stock that had actual investment (not just slightly non-zero values)
    best_stocks(1,count) = varargin{k};
    count = count+1;
end

end

for(k=1:length(best_stocks))

    fprintf("%s\n",best_stocks(k))

end

plot(risk,expected_return)

hold on

plot(risk(lowest_i),expected_return(lowest_i),'o')
end
title("Expected Performance of Stock Portfolio")
xlabel("Risk")
ylabel("Expected Return")
legend("100% limit performance","100% limit alpha","50% limit performance","50% limit alpha","20% limit performance","20% limit alpha","Location","southeast")
axis([0,1.8*(10^-4),0,4.5*(10^-3)])