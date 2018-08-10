close all; clear all; clc;

load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps');
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 's', 'e');
contrastPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/contrastImages/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/';
imageFormat = '*.tif';
imageFiles = dir([contrastPath imageFormat]);
contrastResult = '/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/Curve Fitting results/contrast/';


[rows, cols, ~] = size(AllMaps);
startNum = s + 10;
endNum = e;

for iNum = startNum:endNum
    image = imread([contrastPath imageFiles(iNum).name]);
    map = AllMaps(:,:,iNum - s + 1);
    [x, y] = find(map);
    Poly = createFit(y, x);
    curveX = Poly(y);
    
    figure(1); imshow(image); hold on;
    plot(y, curveX, 'Color', 'y');
    f = getframe;
    imwrite(f.cdata, [contrastResult imageFiles(iNum).name]);
    pause(0.1);
    hold off; 
end

