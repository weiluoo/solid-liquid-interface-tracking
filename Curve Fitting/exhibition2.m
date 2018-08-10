close all; clear all; clc;

originPath = 'F:\DATA-APS\Lianyi_Jul8_2017\80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1\';
imageFormat = '*.tif';
imageFiles = dir([originPath imageFormat]);

s = 97;

% vidObj = VideoWriter('sample.avi');
% vidObj.FrameRate = 10;
% vidObj.Quality = 100;
% open(vidObj);
% rect = [243.51, 113.51, 1501.98, 691.98];

res = {};
startNum = 169;
endNum = 247;
count = 1;
for iNum = startNum:3:endNum
    image = imread([originPath imageFiles(iNum).name]);
    
    % % % % read [ix, iy], [jx, jy] % % %     
    %     [ix, iy] = ?  .mat-(load(sprintf('%d.mat', iNum))), .cvs-(csvread(sprintf('%d.csv', iNum), C1, C2))
    %     [jx, jy] = ?
    iM = csvread(sprintf('%d.csv', iNum));
    iy = iM(:,6); ix = iM(:,7);
    jM = csvread(sprintf('%d.csv', iNum + 3));
    jy = jM(:,6); jx = jM(:,7);
    % % % % % % % % % % % % % 
    
    iPoly = createFit(iy, ix);
    jPoly = createFit(jy, jx);
    iVals = coeffvalues(iPoly);
    iDerivative = iVals(1:9) .* [9:-1:1];
    
    % % % common column % % % 
    commonY = [max(min(iy), min(jy)):0.5:min(max(iy), max(jy))]';
    iXs = iPoly(commonY);
    jXs = jPoly(commonY);
    N = length(commonY);
    % % % % % % % % % % % % % 
    
    % % % compute the normal of line i % % % 
    iValsPoly8 = Poly8(iDerivative, commonY);
    ias = -1 ./ iValsPoly8;
    ibs = iXs - ias .* commonY;
    iNormalVals = repmat(ias, [1, N]) .* repmat(commonY', [N, 1]) + repmat(ibs, [1, N]);
    crossDists = repmat(jXs', [N, 1]) - iNormalVals;
    overZeros = crossDists >= 0;
    % % % % % % % % % % % % % % % % % % % % % 
    
    % % % find the cross points % % % 
    valids = ones(N, 1);
    tmp = sum(overZeros, 2);
    valids = valids - (tmp == 0 | tmp == N);
    [~, idx] = min(abs(crossDists), [], 2);
    crossPoints = [jXs(idx), commonY(idx)];
    % figure(1); hold on;
    % plot(crossPoints(:,2), crossPoints(:,1), 'gx');
    % % % % % % % % % % % % % % % % % 
    
    % % % compute the distance % % %
    dists = sqrt(((crossPoints(:,1) - iXs) .^ 2) + ((crossPoints(:,2) - commonY) .^ 2));
    dists = dists .* valids;
    directs = (crossPoints(:,1) - iXs) > 0;
    % % % % % % % % % % % % % % % % % 
    
    iRes = [iXs, commonY, dists, ias];
    exclude = find(valids == 0);
    iRes(exclude,:) = [];
    res{count} = iRes;
    count = count + 1;
%     csvwrite('Distance.dat', res);
%     type Distance.dat
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1, 2, 1); imshow(image); title(['curve in ', num2str(iNum)]); hold on;
    plot(commonY, iXs, 'Color', 'b', 'LineWidth', 0.1);
    for i = 1:10:N
        if valids(i) == 0
            continue;
        end

        if directs(i) == 1
            plot([commonY(i); crossPoints(i, 2)], [iXs(i); crossPoints(i, 1)], 'Color', 'r', 'LineWidth', 0.1);
        else
            plot([commonY(i); crossPoints(i, 2)], [iXs(i); crossPoints(i, 1)], 'Color', 'g', 'LineWidth', 0.1)
        end
    end
    subplot(1, 2, 2); imshow(image); title(num2str(iNum));
    f = getframe(h);
    imwrite(f.cdata, sprintf('results\\%d.bmp', iNum));
%     writeVideo(vidObj, imcrop(f.cdata, rect));
    pause(0.5);
    close all;
end
save('results.mat', 'res');
% close(vidObj);