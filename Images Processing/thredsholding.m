close all; clear all; clc;

% load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/80_contrast_640.mat');
% load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/161_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1/161_contrast_640.mat');
load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/179_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1/179_contrast_640.mat');

% % % % 80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1 % % % 
% thred = -20;
% limit = 20;
% bound = 8;
% % % % 80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1 % % % 

% % % % 161_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1 % % % 
% thred = -40;
% limit = 20;
% bound = 1;
% % % % 161_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1 % % % 

% % % 179_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1 % % % 
thred = -20;
limit = 10;
bound = 1;
% % % 179_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1 % % % 

[rows, cols, N] = size(F);

for i = 1:N
    f = F(:,:,i);
    BW = f < thred;
    BW1 = bwareaopen(BW, limit);
    D = bwdist(BW1, 'euclidean');
    idx = find(D < bound);
    f(idx) = 0;
%     figure(1); imshowpair(f, F(:,:,i),'montage') 
%     pause(0.1);
    F(:,:,i) = f;
end

% save('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/contrast_thred.mat', 'F', 's', 'e');
% save('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/161_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1/contrast_thred.mat', 'F', 's', 'e');
% save('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries Segment App/Data/179_Al-5Ti-B_t0.4mm_P70S0.5_U18G12_45kHz_100ps_S1/contrast_thred.mat', 'F', 's', 'e');