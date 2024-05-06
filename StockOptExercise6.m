clc
clear variables
close all

start_date = '01012023';
end_date = '02012024';
varargin = {"AMD","KR","TAP","MLM","CXT","COKE","MRO","CMI","WEN","MCD","MKC","PFE","HON","F","JNJ","AMZN","GOOGL","INTC","TXN","PSO","JPM","GS","BAC","LLY","WMT","DIS","PG","ORCL","CMG","BA","AAPL","LUV","HD","GE","GOOG","PEP","CAT","TT","IBM","RTX","LMT","NVDA","MSFT","NFLX","VZ"};
data_raw = hist_stock_data(start_date,end_date,varargin);

n = length(cell2mat(data_raw(3,1,1)));

data = zeros(length(varargin),n);

steps = length(varargin);
tol = 10^-3;
for(i=1:steps)
 P(:,i) = cell2mat(data_raw(3,:,i));

[mu1(i), sigma1(i)] = genMuSigma(P(:,i));

end

%% Calculate Measures

R = (P(2:end,:)-P(1:end-1,:))./P(1:end-1,:);
r = mean(R)';
C = cov(R);

n = length(C);

alpha = [linspace(0,.5,5),1];

selected_stock_indicators = zeros(length(varargin),length(alpha));

for(k=1:length(alpha))
[w(:,k),optval] = quadprog((1-alpha(k)).*2.*C,-alpha(k).*r,[],[],ones(1,n),[1],zeros(n,1),ones(n,1));

expected_return(k) = r'*w(:,k);
selected_stock_indicators(1:length(find(w(:,k)>tol)),k) = find(w(:,k)>tol);



end

for(k=1:length(alpha))
    clear best_stocks
    count = 1;
for(i=1:length(varargin))
if(w(i,k)>=mean(w(:,k)))%find any stock that had actual investment (not just slightly non-zero values)
    best_stocks(1,count) = varargin{i};
    count = count+1;
end
end
fprintf("----alpha = %4.3f----\n",alpha(k))

for(m=1:length(best_stocks))

    fprintf("%s\n",best_stocks(m))

end

end


for(month=1:6)%simulate 4 months and track performance of portfolio, making sure to only take last year of data
    for(j=1:length(varargin))
        Pnew(:,j) = generateStock12_2(P(:,j),20,0);%simulate the next 30 days (abt a month)
        P(:,j) = Pnew(21:end,j);%remove first month of P and add last month of Pnew (the simulated one)
    end

    R = (P(2:end,:)-P(1:end-1,:))./P(1:end-1,:);
    r = mean(R)';
    C = cov(R);
    n = length(C);

    for(k=1:length(alpha))
        expected_return(month,k) = r'*w(:,k);
        [w(:,k),optval] = quadprog((1-alpha(k)).*2.*C,-alpha(k).*r,[],[],ones(1,n),[1],zeros(n,1),ones(n,1));

    end

end

figure

for(k=1:length(alpha))
    
    plot(1:month,expected_return(:,k))
    hold on

end

legend(num2str(alpha'))
xlabel("Months Since Start")
ylabel("Returns")
title("Simulated Returns","2024")
axis([1,month,-1*10^-3,8*10^-3])


%% FUNCTIONS
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
end
end

function [mu, sigma] = genMuSigma(P)

dP = (P(2:end)-P(1:end-1))./P(1:end-1);
mu = mean(dP);
sigma = std(dP);

end
