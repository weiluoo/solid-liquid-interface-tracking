close all; clear all; clc;
% load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/80_contrast_640.mat', 'F', 's', 'e');
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 'F', 's', 'e');
[~, ~, N] = size(F);

for i = 1:N
    figure(1); imshow(F(:,:,i), []);
    f = getframe;
    imwrite(f.cdata, ['threshold/' num2str(s+i-1) '.jpg']);
    pause(0.1);
end

for i = 1:N
    figure(1); imshow(A(:,:,i), []);
    f = getframe;
    imwrite(f.cdata, ['background/' num2str(s+i-1) '.jpg']);
    pause(0.1);
end
