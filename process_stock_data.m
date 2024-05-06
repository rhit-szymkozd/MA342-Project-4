function results = process_stock_data(dataStruct, plt)
%PROCESS_STOCK_DATA Process data from hist_stock_data
%
% inputs:
%         dataStruct: the struct returned from hist_stock_data
%                plt: whether to plot the results
% outputs: 
%            results: 2D matrix

stocks = struct2cell(dataStruct);
numStocks = size(stocks,3);
numVals = 0;

for i = 1:numStocks
    x = stocks(:,:,i);
    x = x(3);
    x = x{1,1};
    j = size(x,1);
    if j > numVals
        numVals = j;
    end
end

results = zeros(numVals, numStocks);

for i = 1:numStocks
    stkVals = stocks(:,:,i);
    vals = stkVals(3);
    vals = vals{1,1};
    results(:,i) = vals;
end

if (plt == 1)
    plot(results);
end