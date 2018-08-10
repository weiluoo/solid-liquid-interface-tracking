close all; clear all; clc;

load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 'F', 's', 'e');
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps');

[rows, cols, n] = size(F);
startNum = 117 - s + 1;
endNum = 125 - s + 1;

for iNum = startNum:endNum
    f = F(:,:,iNum);
    map = AllMaps(:,:,iNum);
    [x, y] = find(map);
    fig = figure(1); imshow(f, []); hold on;
    plot(y, x, 'Color', 'y');
    rect = getrect(fig); rect = round(rect);
    rectangle('Position', rect, 'EdgeColor', 'r');
    
    % % % get the start & end position % % % 
    startP = rect(1) - y(1) + 1;
    endP = length(y);
    lenX = 2; shiftX = (-lenX:0);
    
    % % % get the costMap % % % 
    costMap = -f;
    costMap = costMap .* (costMap > 0);
    h = fspecial('gaussian', 5, 1);
    costMap = xcorr2(costMap, h);
    costMap = costMap(3:end-2,:);
    costMap = costMap(:,3:end-2);
    
    % % % initialize % % % 
    paths = x(1:startP-1);
    thred = 4;
    currCost = -inf;
    
    for p = startP:endP
        col = y(p);
        prevRow = paths(p-1);
        currRows = prevRow + shiftX;
        [~, idxs] = sort(costMap(currRows ,col), 'descend');
        currRows = currRows(idxs);
        curvatures = (mean(paths(p-7:p-5)) - 2 * mean(paths(p-4:p-2)) + currRows) .^ 2;
        newCost = costMap(currRows, col) * (1 - curvatures ./ thred);
        for i = 1:3
            if curvatures(i) < thred & newCost(i) > currCost
                paths = [paths; currRows(i)];
                currCost = newCost(i);
                break;
            end
        end
    end

    plot(y, paths, 'Color', 'r');
    pause(0.1);
    
    newMap = zeros(rows, cols);
    for i = 1:length(y)
        newMap(paths(i), y(i)) = 1;
    end
    AllMaps(:,:,iNum) = newMap;
end

% save('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps');