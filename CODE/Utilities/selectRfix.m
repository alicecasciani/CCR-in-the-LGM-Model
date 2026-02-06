function R = selectRfix(modelType, S_black, S_bach)
% selectRfix Returns fixed rate vector based on model type
%   R = selectRfix(modelType, S_black, S_bach) picks swap rates
%   corresponding to the chosen volatility model.
%
%   INPUTS:
%     modelType - 'BLACK' or 'BACHELIER'
%     S_black   - Vector of swap rates for Black model
%     S_bach    - Vector of swap rates for Bachelier model
%
%   OUTPUT:
%     R         - Selected swap rate vector
    switch modelType
        case 'BLACK'
            R = S_black;
        case 'BACHELIER'
            R = S_bach;
    end
end