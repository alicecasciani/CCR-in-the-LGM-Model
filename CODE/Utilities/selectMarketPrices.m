function prices = selectMarketPrices(modelType, blackPrices, bachPrices)
% selectMarketPrices Returns market price vector based on model type
%   prices = selectMarketPrices(modelType, blackPrices, bachPrices)
%   chooses between Black and Bachelier price arrays.
%
%   INPUTS:
%     modelType   - 'BLACK' or 'BACHELIER'
%     blackPrices - Vector of observed Black-model prices
%     bachPrices  - Vector of observed Bachelier-model prices
%
%   OUTPUT:
%     prices      - Selected price vector
    switch modelType
        case 'BLACK'
            prices = blackPrices;
        case 'BACHELIER'
            prices = bachPrices;
        otherwise
            error('Unknown model type: %s', modelType);
    end
end