close all; clear all; clc;

load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 'F', 's', 'e');
% load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/161_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1/contrast_thred.mat', 'F', 's', 'e');
% load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/179_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1/contrast_thred.mat', 'F', 's', 'e');


originPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/contrastImages/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/';
% originPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/contrastImages/161_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1/';
% originPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/contrastImages/179_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1/';
imageFormat = '*.tif';
imageFiles = dir([originPath imageFormat]);


index = 80;
image = imread([originPath imageFiles(s + index - 1).name]);
f = F(:,:,index);
[rows, cols] = size(f);


h = fspecial('gaussian', 5, 1);
cost = -f .* (f < 0);
cost = xcorr2(cost, h);
cost = cost(3:end-2,:);
cost = cost(:,3:end-2);


fig = figure(1); imshow(image); 
% [x, y] = getline(fig);
[y, x] = getpts(fig);
y = round(y); x = round(x);
hold on;
% plot(y, x, 'yx');


map = zeros(rows, cols);
N = length(y);
y1 = y(1); x1 = x(1);
plot(y1, x1, 'yx', 'markers', 12);
for i = 2:N
    y2 = y(i); x2 = x(i);
    plot(y2, x2, 'yx', 'markers', 12);
    minRow = min(x1, x2) - 8;
    maxRow = max(x1, x2) + 8;
    costMap = cost(minRow:maxRow, y1:y2);
    [shortestPath, pathMap] = DP(costMap, x1 - minRow + 1, x2 - minRow + 1);
    map(minRow:maxRow, y1:y2) = pathMap;
    plot(shortestPath(:,2) + y1 - 1, shortestPath(:,1) + minRow - 1, 'Color', 'w');
    y1 = y2; x1 = x2;
end


% save('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/demo.mat', 'map', 'index', 's');
