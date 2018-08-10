function [shortestPath, pathMap] = DP(costMap, startRow, endRow)
    lenX = 2; shiftX = (-lenX:lenX)'; weights = exp(-(-lenX:lenX) .^2 / (5 * 25))';
    [nrows, ncols] = size(costMap);
    
    C = ones(nrows, ncols) .* inf;
    C(startRow, 1) = 0;
    T = zeros(nrows, ncols);
    
    for col = 2:ncols
        for row = 1:nrows
            prevRows = row + shiftX;
            w = weights;
            
            exclude = (prevRows < 1) | (prevRows > nrows);
            prevRows(exclude) = [];
            w(exclude) = [];
            exclude = C(prevRows, col - 1) == inf;
            prevRows(exclude) = [];
            w(exclude) = [];
            
            if isempty(prevRows)
                continue;
            end
            [minC, idx] = max(C(prevRows, col - 1) .* w);
            C(row, col) = minC + costMap(row, col);
            T(row, col) = prevRows(idx);
        end
    end
    
    pathMap = zeros(nrows, ncols);
    row = endRow;
    pathMap(row, ncols) = 1;
    shortestPath = [row, ncols];
    for col = ncols:-1:2
        row = T(row, col);
        pathMap(row, col - 1) = 1;
        shortestPath = [row, col - 1 ;shortestPath];
    end
end