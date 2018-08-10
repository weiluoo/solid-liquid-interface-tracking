close all; clear all; clc;

load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 'F', 's', 'e');
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps', 'index');

[rows, cols, n] = size(F);
bound = 390;
startColMove = 10;
endColMove = 8;
for iNum = index+1:n
    f = F(:,:,iNum);
    prevMap = AllMaps(:,:,iNum-1);
    
    % % % % % dynamic programming % % % % % 
    thred = 10;
    
    % % % get the mask % % % 
    eudist = bwdist(prevMap, 'euclidean');
    mask = eudist < thred;
    mask(1:bound-1,:) = 0;
    
    % % get the cost map % % % 
    costMap = -f .* mask;
    costMap = costMap .* (costMap > 0);
    h = fspecial('gaussian', 5, 1);
    costMap = xcorr2(costMap, h);
    costMap = costMap(3:end-2,:);
    costMap = costMap(:,3:end-2);
    
    % % % find the startCol and endCol % % %     
    [prevX, prevY] = find(prevMap);
    startArea = costMap(prevX(1)-1:prevX(1), prevY(1)+1:prevY(1)+max(startColMove+2, 10));
    [~, idxY] = max(sum(startArea, 1));
    startCol = prevY(1) + idxY;
    startColMove = idxY;
    if prevY(end)+10 >= cols
        endCol = cols;
    else
        endArea = costMap(prevX(end)-5:prevX(end), prevY(end)+1:prevY(end)+max(endColMove+2, 10));
        [~, idxY] = max(sum(endArea, 1));
        endCol = prevY(end) + idxY;
        endColMove = idxY;
    end
    
    minRow = min(prevX) - thred; maxRow = max(prevX) + thred;
    costMap = costMap(minRow:maxRow, startCol:endCol);
    mask = mask(minRow:maxRow, startCol:endCol);
    
    % % % initalize % % %
    [nrows, ncols] = size(mask);
    paths = [];
    
    % % % for the start col % % %
    for row = 1:nrows
        if mask(row, 1) == 0
            continue;
        end
        paths = [paths; costMap(row, 1), row];
    end
    
    % % % for the second col % % %
    lenX = 2; shiftX = (-lenX:lenX)';
    newPaths = [];
    for row = 1:nrows
        if mask(row, 2) == 0
            continue;
        end
        prevRows = row + shiftX;
        % remove the points outside the picture
        exclude = (prevRows < 1) | (prevRows > nrows);
        prevRows(exclude) = [];
        % remove the points which are inactivate
        exclude = mask(prevRows, 2) == 0;
        prevRows(exclude) = [];
        if isempty(prevRows)
            mask(row, 2) = 0;
            continue;
        end
        pathsIdx = findIdxs(paths(:,end), prevRows);
        availablePaths = paths(pathsIdx,:);
        availablePaths = [availablePaths, ones(size(availablePaths, 1), 1) * row];
        availablePaths(:,1) = availablePaths(:,1) + costMap(row, 2);
        newPaths = [newPaths; availablePaths];
    end
    paths = newPaths;
    
    % % % for the rest cols % % %
    for col = 3:ncols
        newPaths = [];
        for row = 1:nrows
            if mask(row, col) == 0
                continue;
            end
            prevRows = row + shiftX;
            % remove the points outside the picture
            exclude = (prevRows < 1) | (prevRows > nrows);
            prevRows(exclude) = [];
            % remove the points which are inactivate
            exclude = mask(prevRows, col) == 0;
            prevRows(exclude) = [];
            if isempty(prevRows)
                mask(row, col) = 0;
                continue;
            end
            
            for j = 1:size(prevRows, 1)
                pathsIdx = findIdxs(paths(:,end), prevRows(j));
                availablePaths = paths(pathsIdx,:);
                availablePaths = [availablePaths, ones(size(availablePaths, 1), 1) * row];
                % % % curvature % % %
                curvatures = (availablePaths(:,end-2) - 2 * availablePaths(:,end-1) + availablePaths(:,end)) .^ 2;
                curvatures = sqrt(curvatures) / 4;
                curvatures = max(curvatures, (curvatures >= 0.5) * 100) / 2;
                % % % % % % % % % % %
                availablePaths(:,1) = availablePaths(:,1) + costMap(row, col) - curvatures;
                [maxVal maxIdx] = max(availablePaths(:,1));
                newPaths = [newPaths; availablePaths(maxIdx,:)];
            end
        end
        paths = newPaths;
    end
    
    % % % find the highest scores path % % % 
    [~, Idx] = max(paths(:,1));
    figure(1); imshow(f, []); title(num2str(iNum)); hold on;
    plot(startCol:endCol, paths(Idx, 2:end)+minRow, 'Color', 'y');
    hold off;
    
    iMap = zeros(rows, cols);
    for col = startCol:endCol
        iMap(paths(Idx, col-startCol+2)+minRow, col) = 1;
    end

    AllMaps(:,:,iNum) = iMap;
    pause(0.1);
end

% save('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps', 'index', 'bound');