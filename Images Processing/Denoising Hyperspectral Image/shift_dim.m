clc; %clear the command window 
close all; %close all figure windows
clear all; %clear all variables in the workspace

Hyp_dark = multibandread('SBDay2Dark_VNIR.raw',[100,320,256],'float',0,'bip','ieee-le');
Hyp_image = multibandread('SBDay2Scan2_VNIR.raw',[3868,320,256],'float',0,'bip','ieee-le');
Hyp_dark_mean = mean(Hyp_dark);
Hyp = Hyp_image-repmat(Hyp_dark_mean,[3868,1]);

I = Hyp(:,:,125);
J = shiftdim(I,1);
imshow(J,[]);
