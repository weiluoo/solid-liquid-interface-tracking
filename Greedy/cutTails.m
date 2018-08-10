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
    minCol = rect(1);
    
    map(:,minCol:end) = 0;
    [x, y] = find(map);
    plot(y, x, 'Color', 'r');
    
    AllMaps(:,:,iNum) = map;
end

% save('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps');