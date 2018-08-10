close all; clear all; clc;
warning off;

load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps');
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 's', 'e');
contrastPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/contrastImages/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/';
imageFormat = '*.tif';
imageFiles = dir([contrastPath imageFormat]);

groundtruePath = 'groundTrue/';
startNum = 178;
endNum = 238;

e = [];

for iNum = startNum:3:endNum
    image = imread([contrastPath imageFiles(iNum).name]);
    map = AllMaps(:,:,iNum - s + 2);
    [x, y] = find(map);
    Poly = createFit(y, x);
    minY = y(1); maxY = y(end);
    gt = csvread([groundtruePath num2str(iNum) '.csv']);
    x = gt(:,7); y = gt(:,6);
    gtPoly = createFit(y, x);
    minY = max(minY, round(y(1))) + 5;
    maxY = min(maxY, round(y(end))) -5;
    y = [minY:0.1:maxY]';
    PolyX = Poly(y);
    gtPolyX = gtPoly(y);
%     figure(1); imshow(image); hold on;
%     plot(y, gtPolyX, 'Color', 'r', 'LineWidth', 1);
%     f = getframe;
%     imwrite(f.cdata, ['example/gt_' num2str(iNum) '.jpg']);
%     hold off;
%     figure(2); imshow(image); hold on;
%     plot(y, PolyX, 'Color', 'y', 'LineWidth', 1);
%     f = getframe;
%     imwrite(f.cdata, ['example/' num2str(iNum) '.jpg']);
%     hold off;
    ie = mean((abs(PolyX - gtPolyX) ./ gtPolyX) * 100);
    e = [e; ie];
    disp(num2str(ie));
    pause(0.1);
end

% 0.4202
% 0.31435
% 0.70072
% 0.43525
% 0.55784
% 0.78258
% 0.65787
% 1.2041
% 0.3807
% 0.51916
% 0.5984
% 0.36465
% 1.8156
% 0.17837
% 3.6867
% 0.44057
% 1.1474
% 0.78487
% 0.65933
% 0.57343
% 1.2362