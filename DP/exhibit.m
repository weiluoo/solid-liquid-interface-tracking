close all; clear all; clc;

load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/AllMaps.mat', 'AllMaps');
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 'F', 's', 'e');
FResult = '/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/Greedy results/F/';

originPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/origin/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/';
imageFormat = '*.tif';
originFiles = dir([originPath imageFormat]);
originResult = '/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/Greedy results/origin/';

contrastPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/contrastImages/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/';
contrastFiles = dir([contrastPath imageFormat]);
contrastResult = '/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Results/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/Greedy results/contrast/';

for iNum = s+9:e
    
    originImage = imread([originPath originFiles(iNum).name]);
    contrastImage = imread([contrastPath contrastFiles(iNum).name]);
    i = iNum - s + 1;
    f = F(:,:,i);
    map = AllMaps(:,:,i);
    [x, y] = find(map);
    
    h1 = figure(1); imshow(originImage); hold on;
    plot(y, x, 'Color', 'y');
    f1 = getframe(h1);
    image1 = f1.cdata; imwrite(image1(29:668, 93:716,:), [originResult originFiles(iNum).name]);
    hold off;
    
    h2 = figure(2); imshow(contrastImage); hold on;
    plot(y, x, 'Color', 'y');
    f2 = getframe(h2);
    image2 = f2.cdata; imwrite(image2(29:668, 93:716,:), [contrastResult contrastFiles(iNum).name]);
    hold off;
    
    h3 = figure(3); imshow(f, []); hold on;
    plot(y, x, 'Color', 'y');
    f3 = getframe(h3);
    image3 = f3.cdata; imwrite(image3(29:668, 93:716,:), [FResult contrastFiles(iNum).name]);
    hold off;
    
    pause(0.2);
end