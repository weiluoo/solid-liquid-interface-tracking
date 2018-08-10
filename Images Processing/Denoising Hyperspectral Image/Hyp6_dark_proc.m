%
%load dak scan ('file_name',[lines,samples,bands],'uint16',0,'bil','ieee-le');
Hyp_dark_1=multibandread('SBDay2Dark_SWIR.raw',[100,320,256],'uint16',0,'bil','ieee-le');
Hyp_dark_2=multibandread('SBDay2Dark_VNIR.raw',[100,320,256],'uint16',0,'bil','ieee-le');
Hyp_dark_3=multibandread('SBDay2Scan2_VNIR.raw',[100,320,256],'uint16',0,'bil','ieee-le');
Hyp_dark_4=multibandread('SBDay2Scan3_SWIR.raw',[100,320,256],'uint16',0,'bil','ieee-le');
for i=1:256
m=Hyp_dark_4(:,:,i);
figure(i);
imshow(m,[]);
end

%%
%load dak scan ('file_name',[lines,samples,bands],'uint16',0,'bil','ieee-le');
%Hyp_dark=multibandread('SBDay2Dark_SWIR.raw',[100,320,256],'uint16',0,'bil','ieee-le');
%%
% %average the line. The '100' below may need to be changed.
% Hyp_dark_mean=(sum(Hyp_dark))/100;
% %% Use different approach
% Hyp_dark_mean1=mean(Hyp_dark);
% %%
% %Import the image to be corrected, again, edit the number of lines.
% Hyp_Image=multibandread('SBDay2Scan3_SWIR.raw',[3868,320,256],'uint16',0,'bil','ieee-le');
% %%
% %subtract the averaged dark scan scene. 
% %‘W’ below allows depicting the corrected image. 
% Destripe_Data=bsxfun(@minus,Hyp_Image,Hyp_dark_mean);
% Destripe_Data1=bsxfun(@minus,Hyp_Image,Hyp_dark_mean1);
% %W = squeeze(Img1_drcr(1:50,1:320,50));
% %imagesc(W)
% %%
% %Export the destriped image to a new .raw file
% multibandwrite(Destripe_Data,'SBDay2Scan3_SWIR_destrip.raw','bil');
% multibandwrite(Destripe_Data1,'SBDay2Scan3_SWIR_destrip1.raw','bil');
% %%
