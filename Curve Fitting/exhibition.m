close all; clear all; clc;
warning off;
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps');
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 's', 'e');
originPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/contrastImages/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/';
imageFormat = '*.tif';
imageFiles = dir([originPath imageFormat]);
contrastResult = '/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/Curve Fitting results/origin/';


res = {};
count = 1;
% h = figure('units','normalized','outerposition',[0 0 1 1]);
for iNum = s+15:e-10
    image = imread([originPath imageFiles(iNum).name]);
    iTrace = AllMaps(:,:,iNum-s+1);
    jNum = iNum + 2;
    jTrace = AllMaps(:,:,jNum-s+1);
    [ix, iy] = find(iTrace);
    [jx, jy] = find(jTrace);
    
    iPoly = createFit(iy, ix);
    jPoly = createFit(jy, jx);
    iVals = coeffvalues(iPoly);
    iDerivative = iVals(1:9) .* [9:-1:1];
    
    % % % common column % % % 
    commonY = [max(min(iy), min(jy)):0.1:min(max(iy), max(jy))]';
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
    
    figure(1);
%     h = figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(image); hold on;
    plot(commonY, iXs, 'Color', 'w', 'LineWidth', 0.1);
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
%     f = getframe;
%     imwrite(f.cdata, [contrastResult imageFiles(iNum).name]);
%     writeVideo(vidObj, imcrop(f.cdata, rect));
%     pause(0.1);
%     close all;
    hold off;
    pause(0.1);
end

% close(vidObj);