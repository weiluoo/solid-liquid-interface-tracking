close all; clear all; clc;

load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 'F', 's', 'e');
% load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps');
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/DP results/AllMaps.mat', 'AllMaps');

contrastPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/contrastImages/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/';
imageFormat = '*.tif';
contrastFiles = dir([contrastPath imageFormat]);

iNum = 189;
image = imread([contrastPath contrastFiles(iNum).name]);
iNum = iNum - s + 1;
count = 0;

f = F(:,:,iNum);
map = AllMaps(:,:,iNum);
[x, y] = find(map);
[rows, cols] = size(f);

whiteboard = ones(rows, cols);

fig = figure(1); imshow(image, []); hold on;
plot(y, x, 'Color', 'y', 'LineWidth', 1);
% rect = getrect(fig); rect = round(rect);
rect0 = [117, 388, 423, 138];
rect = [408, 468, 77, 36];
rectangle('Position', rect, 'EdgeColor', 'k');
ff = getframe;
imwrite(imcrop(ff.cdata, rect0), ['examples/', num2str(count), '.jpg']);

% % % get the start & end position % % % 
startP = rect(1) - y(1) + 1;
endP = rect(1) + rect(3) - y(1) + 1;
minRow = rect(2); maxRow = rect(2) + rect(4);
lenX = 2; shiftX = (-lenX:lenX);

% % % get the costMap % % % 
mask = zeros(rows, cols);
mask(minRow:maxRow, rect(1):rect(1) + rect(3)) = 1;
costMap = -f .* mask;
costMap = costMap .* (costMap > 0);
h = fspecial('gaussian', 5, 1);
costMap = xcorr2(costMap, h);
costMap = costMap(3:end-2,:);
costMap = costMap(:,3:end-2);

% % % initialize % % % 
paths = [x(1:startP-1)', pathGenerator(x(startP), x(endP), endP-startP), x(endP+1:end)']';
curvatureS = (mean(paths(startP-7:startP-5)) - 2 * mean(paths(startP-4:startP-2)) + paths(startP)) .^ 2;
curvatureE = (mean(paths(endP+5:endP+7)) - 2 * mean(paths(endP+2:endP+4)) + paths(endP)) .^ 2;
d = y(endP) - y(startP) - (abs(paths(endP) - paths(startP)) / 2);
minRow = max(minRow, min(paths(startP), paths(endP)) - d);
maxRow = min(maxRow, max(paths(startP), paths(endP)) + d);
thred = 3.5;

while endP > startP
    count = count + 1;
    hold off; imshow(whiteboard, []); hold on;
    rectangle('Position', rect, 'EdgeColor', 'k');
    plot(y, paths, 'Color', 'b', 'LineWidth', 1);
    plot(y(startP:endP), paths(startP:endP), 'Color', 'r', 'LineWidth', 1);
    ff = getframe;
    imwrite(imcrop(ff.cdata, rect0), ['examples/' num2str(count) '.jpg']);
    if curvatureS >= curvatureE
        prevRow = paths(startP-1);
        currRows = prevRow + shiftX;
        % % % out of search area % % %
        exclude = (currRows < minRow) | (currRows > maxRow);
        currRows(exclude) = [];
        % % % curvatures for previous point % % %
        prevCurvatures = (paths(startP-2) - 2 * paths(startP-1) + currRows) .^ 2;
        prevCurvatures = sqrt(prevCurvatures) / 4;
        exclude = (prevCurvatures > 0.5);
        currRows(exclude) = [];
        % % % curvatures for the endP and startP % % %
        [~, idxs] = sort(costMap(currRows ,y(startP)), 'descend');
        currRows = currRows(idxs);
        d = std(costMap(currRows ,y(startP)));
        curvatureE = (mean(paths(endP-4:endP-2)) - 2 * paths(endP) + mean(paths(endP+2:endP+4))) .^ 2;
        currCost = -inf;
        for i = 1:size(currRows, 2)
            row = currRows(i);
            newPath = [paths(1:startP-1)', pathGenerator(row, paths(endP), endP-startP), paths(endP+1:end)']';
            newS = (mean(newPath(startP-7:startP-5)) - 2 * mean(newPath(startP-4:startP-2)) + newPath(startP)) .^ 2;
            newE = (mean(newPath(endP-4:endP-2)) - 2 * newPath(endP) + mean(newPath(endP+2:endP+4))) .^ 2;
            newCost = (1 - newS / thred) * costMap(row, y(startP));
            if (newS < thred) & (newE < thred) & (newCost > currCost)
                paths = newPath;
                curvatureS = newS;
                curvatureE = newE;
                currCost = newCost;
                break;
            end
        end
        startP = startP + 1;
    else
        prevRow = paths(endP+1);
        currRows = prevRow + shiftX;
        % % % out of search area % % %
        exclude = (currRows < minRow) | (currRows > maxRow);
        currRows(exclude) = [];
        % % % curvatures for previous point % % %
        prevCurvatures = (paths(endP+2) - 2 * paths(endP+1) + currRows) .^ 2;
        prevCurvatures = sqrt(prevCurvatures) / 4;
        exclude = (prevCurvatures > 0.5);
        currRows(exclude) = [];
        % % % curvatures for the endP and startP % % %
        [~, idxs] = sort(costMap(currRows ,y(endP)), 'descend');
        currRows = currRows(idxs);
        d = std(costMap(currRows ,y(endP)));
        curvatureS = (mean(paths(startP-4:startP-2)) - 2 * paths(startP) + mean(paths(startP+2:startP+4))) .^ 2;
        currCost = -inf;
        for i = 1:size(currRows, 2)
            row = currRows(i);
            newPath = [paths(1:startP-1)', pathGenerator(paths(startP), row, endP-startP), paths(endP+1:end)']';
            newE = (mean(newPath(endP+5:endP+7)) - 2 * mean(newPath(endP+2:endP+4)) + newPath(endP)) .^ 2;
            newS = (mean(newPath(startP-4:startP-2)) - 2 * newPath(startP) + mean(newPath(startP+2:startP+4))) .^ 2;
            newCost = (1 - newE / thred) * costMap(row, y(endP));
            if (newS < thred) & (newE < thred) & (newCost > currCost)
                paths = newPath;
                curvatureE = newE;
                curvatureS = newS;
                currCost = newCost;
                break;
            end
        end
        endP = endP - 1;
    end
    curvatureS = (mean(paths(startP-7:startP-5)) - 2 * mean(paths(startP-4:startP-2)) + paths(startP)) .^ 2;
    curvatureE = (mean(paths(endP+5:endP+7)) - 2 * mean(paths(endP+2:endP+4)) + paths(endP)) .^ 2;
    d = y(endP) - y(startP) - (abs(paths(endP) - paths(startP)) / 2);
    minRow = max(minRow, min(paths(startP), paths(endP)) - d);
    maxRow = min(maxRow, max(paths(startP), paths(endP)) + d);
    pause(0.1)
end

newMap = zeros(rows, cols);
for i = 1:length(y)
    newMap(paths(i), y(i)) = 1;
end
AllMaps(:,:,iNum) = newMap;

% save('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps');
