% This is the modified code of the main code. In this code, we add an item
% to the objective function.
clc; %clear the command window 
close all; %close all figure windows
clear all; %clear all variables in the workspace

addpath siims_release;

load('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries/rmbk_316.mat', 'F');

[rows, cols, n] = size(F);
Hyp = F;

% Read the images and take a subimage of 160*624*1249
Hyp_actual = double(Hyp);
[row, col, band] = size(Hyp_actual);
Hyp_actual_normalized = zeros(row, col, band);
Observe = zeros(row*col, band);

% Create a vector to store the maximum value of each band
Max = zeros(1, band);
Min = zeros(1, band);

% Create a vector to store the random Gaussian noise variance.
% Parameter = 0.5*rand(1,band);
for g = 1:band
    Max(g) = max(max(Hyp_actual(:,:,g)));
    Min(g) = min(min(Hyp_actual(:,:,g)));
    Hyp_actual_normalized(:,:,g) = (Hyp_actual(:,:,g) - Min(g))./(Max(g)-Min(g));
end

weight = ones(1, band);


%-------------------------------------------------------------------------------------------------------------------------
% Denoising the noisy images

% Vectorizing each band
for j = 1:band
    m = Hyp_actual_normalized(:,:,j);
    n = m(:);
    Observe(:,j) = n;
end


[A_hat, A_h, ~] = inexact_alm_rpca(Observe); % original one, without TV_norm
[F_hat, F_h, ~] = LRTV_11(Observe, 1/316, 0.01*weight, 1e-8, 1000, row, col); % with TV_norm



% Reshape the restored images
A = zeros(row, col, band);
F = zeros(row, col, band);

% for k = 1:band
%     a = A_hat(:,k);
%     AA = reshape(a,row, col);
%     A(:,:,k) = AA;
%     A(:,:,k) = A(:,:,k).*(Max(k)-Min(k)) + Min(k);
% 
%     f = F_hat(:,k);
%     FF = reshape(f, row, col);
%     F(:,:,k) = FF;
%     F(:,:,k) = F(:,:,k).*(Max(k)-Min(k)) + Min(k);
% end

for k = 1:band
    a = F_hat(:,k);
    AA = reshape(a,row, col);
    A(:,:,k) = AA;
    A(:,:,k) = A(:,:,k).*(Max(k)-Min(k)) + Min(k);

    f = F_h(:,k);
    FF = reshape(f, row, col);
    F(:,:,k) = FF;
    F(:,:,k) = F(:,:,k).*(Max(k)-Min(k)) + Min(k);
end

%% =====================================================================================

% bandth2 = 20;
% 
% figure(1);
% imshow(Hyp_actual(:,:,bandth2),[]);
% 
% figure(2);
% imshow(A(:,:,bandth2),[]);
% 
% figure(3);
% imshow(F(:,:,bandth2),[]);

% vidObj = VideoWriter('currentResult.avi');
% vidObj.FrameRate = 5;
% vidObj.Quality = 99;
% 
% open(vidObj);

% for bandth3 = 1:band
%     figure(1);
%     subplot(3, 1, 1); imshow(Hyp_actual(:,:,bandth3),[]); title('original image');
%     subplot(3, 1, 2); imshow(A(:,:,bandth3),[]); title('background');
%     subplot(3, 1, 3); imshow(F(:,:,bandth3),[]); title('boundaries');
%     pause(0.1);
% %     f = getframe(figure(1));
% %     writeVideo(vidObj, f.cdata);
% end
% close(vidObj);

save('/usr/local/home/wlhz2/Documents/MATLAB/Boundaries/rmbk_316.mat', 'F');

