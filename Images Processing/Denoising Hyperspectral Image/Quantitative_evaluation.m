% This is the modified code of the main code. In this code, we add an item
% to the objective function.
clc; %clear the command window 
close all; %close all figure windows
clear all; %clear all variables in the workspace

addpath siims_release;

% Parameters
sigma = 0.05; %variance of Gaussian noise
den = 0.10; %noise density of Salt and Pepper noise
alpha = 0.01; % tradeoff the weight of the smooth term.
Size = 5; % size of the filter window
% band = 191; % number of bands

% Read the images and take a subimage of 256*256*191
Hyp = imread('dc.tif');
Hyp_actual = Hyp(531:786,1:256,:);
Hyp_actual = double(Hyp_actual);
[row, col, band] = size(Hyp_actual);
Hyp_noise = zeros(row, col, band);
Hyp_actual_normalized = zeros(row, col, band);
Observe = zeros(row*col,band);

% Create a vector to store the maximum value of each band
Max = zeros(1,band);
Min = zeros(1,band);

% Create a vector to store the random Gaussian noise variance.
% Parameter = 0.5*rand(1,band);
Parameter = 0.1;
parameter = zeros(1,band);
for g = 1:band
    Max(g) = max(max(Hyp_actual(:,:,g)));
    Min(g) = min(min(Hyp_actual(:,:,g)));
    Hyp_actual_normalized(:,:,g) = (Hyp_actual(:,:,g) - Min(g))./(Max(g)-Min(g));
end

% Add noise to every band and normalize each band
% H = fspecial('Gaussian', [11 11], sigma); % Gaussian lowpass filter
for i = 1:band
    Hyp_noise(:,:,i) = imnoise(Hyp_actual_normalized(:,:,i), 'gaussian' , 0, sigma); % add the Gaussian noise randomly
%     Hyp_noise(:,:,i) = Hyp_actual_normalized(:,:,i) + Parameter(i)*randn(size(Hyp_actual_normalized(:,:,i))); % add the Gaussian noise randomly
    Hyp_noise(:,:,i) = imnoise(Hyp_noise(:,:,i),'salt & pepper', Parameter); % add the salt and pepper noise randomly
%     Hyp_noise(100,101:200,i) = zeros(1,100); % add some strip noise at the same position of each band
    Hyp_noise(151:152, 51:150,i) = zeros(2,100); % add some strip noise at the same position of each band
end

mask = zeros(row, col);
% mask(100, 101:200) = ones(1,100);
mask(151:152, 51:150) = ones(2,100);

column_number = randi(30, 1, band); % number of bands chosen to add noise
for j = 1:band
    choose_column = randi([1,256],1,column_number(j));
    for l = 1:length(choose_column)
        Hyp_noise(:,choose_column(l),j) = zeros(256,1); % add the stripe noise randomly
    end
end

% Calculate the gradient of each band as the reference of the weights of HTV term
% weight = zeros(1,band);
% for k = 1:band
%     [dx, dy] = gradient(Hyp_noise(:,:,k));
%     weight(k) = sum(sum(sqrt(dx.^2+dy.^2)));
% end
% weight = ((weight - min(weight))/(max(weight) - min(weight)) + 1)/2; % normalize the weights

weight = ones(1, band);

% Hyp_noise(:,:,126:135) = zeros(256,256,10); % corrupt the whole image of band 126-135

% % calculate the parameter, i.e. the similarity between two adjacent bands
% for ii = 2:band
%     parameter(ii) = alpha * SSIM(Hyp_noise(:,:,ii),Hyp_noise(:,:,ii-1),[0.01 0.03],fspecial('gaussian', 11, 1.5),1);
% end
% parameter(1) = alpha;
 
%-------------------------------------------------------------------------------------------------------------------------
% Denoising the noisy images

% Denoise hyperspectral images with filter functions
% X_hat = zeros(256^2,band);
% Y_hat = zeros(256^2,band);
% Z_hat = zeros(256^2,band);
% 
% for jj = 1:band
%     image = Hyp_noise(:,:,jj);
%     h1 = fspecial('average',Size); % average filter function
%     h2 = fspecial('gaussian',Size,0.5); % Gaussian filter function
%     denoise_image_1 = medfilt2(image,[Size,Size]); % denoise with median filter
%     denoise_image_2 = filter2(h1, image); % denoise with average filter
%     denoise_image_3 = filter2(h2, image); % denoise with Gaussian filter
%     X_hat(:,jj) = denoise_image_1(:);
%     Y_hat(:,jj) = denoise_image_2(:);
%     Z_hat(:,jj) = denoise_image_3(:);
% end

% Vectorizing each band
for j = 1:band
    m = Hyp_noise(:,:,j);
    n = m(:);
    Observe(:,j) = n;
end

% % Construct the smooth matrices A_smooth and B_smooth
% ncols = size(Observe,2);
% % A_smooth = zeros(256,256);
% B_smooth = zeros(ncols,ncols);
% % % A_smooth(1,1:2) = [.5,.5]; A_smooth(256,255:256) = [.5,.5];
% % B_smooth(1:2,1) = [.5,.5]; B_smooth(ncols-1:ncols,ncols) = [.5,.5];
% % 
% % % for s = 2:255
% % %     A_smooth(s,s-1:s+1) = [1/3,1/3,1/3];
% % % end
% % 
% % for r = 2:ncols-1
% %     B_smooth(r-1:r+1,r) = [1/3, 1/3, 1/3];
% % end
% 
% for r = 1:ncols
%     B_smooth(r,r) = 1;
% end

% for r = 124:137
%     B_smooth(r-1:r+1,r) = [1/3,1/3,1/3];
% end

%Augmented Lagrange Multiplier (ALM)-inexact


[A_hat, A_h, ~] = inexact_alm_rpca(Observe); % original one, without TV_norm

% [B_hat, ~] = LRTV_11(Observe, 1/256, 0.001*weight, 1e-8, 1000);
% [C_hat, ~] = LRTV_11(Observe, 1/256, 0.002*weight, 1e-8, 1000);
% [D_hat, ~] = LRTV_11(Observe, 1/256, 0.005*weight, 1e-8, 1000);
% [E_hat,~] = LRTV_11(Observe, 1/256, 0.01*weight, 1e-8, 1000);
[F_hat, ~] = LRTV_11(Observe, 1/256, 0.01*weight, 1e-8, 1000);
% [G_hat, ~] = LRTV_11(Observe, 1/256, 0.05*weight, 1e-8, 1000);
% [H_hat, ~] = LRTV_11(Observe, 1/256, 0.1*weight, 1e-8, 1000);


% Reshape the restored images
A = zeros(256,256,band);
% B = zeros(256,256,band);
% C = zeros(256,256,band);
% D = zeros(256,256,band);
% E = zeros(256,256,band);
F = zeros(256,256,band);
% G = zeros(256,256,band);
% H = zeros(256,256,band);
% X = zeros(256,256,band);
% Y = zeros(256,256,band);
% Z = zeros(256,256,band);
for k = 1:band
    a = A_hat(:,k);
    AA = reshape(a,256,256);
    A(:,:,k) = AA;
    A(:,:,k) = A(:,:,k).*(Max(k)-Min(k)) + Min(k).*ones(256,256);
%     b = B_hat(:,k);
%     BB = reshape(b,256,256);
%     B(:,:,k) = BB;
%     B(:,:,k) = B(:,:,k).*(Max(k)-Min(k)) + Min(k).*ones(256,256);
%     c = C_hat(:,k);
%     CC = reshape(c,256,256);
%     C(:,:,k) = CC;
%     C(:,:,k) = C(:,:,k).*(Max(k)-Min(k)) + Min(k).*ones(256,256);
%     d = D_hat(:,k);
%     DD = reshape(d,256,256);
%     D(:,:,k) = DD;
%     D(:,:,k) = D(:,:,k).*(Max(k)-Min(k)) + Min(k);
%     e = E_hat(:,k);
%     EE = reshape(e,256,256);
%     E(:,:,k) = EE;
%     E(:,:,k) = E(:,:,k).*(Max(k)-Min(k)) + Min(k);
    f = F_hat(:,k);
    FF = reshape(f,256,256);
    F(:,:,k) = FF;
    F(:,:,k) = F(:,:,k).*(Max(k)-Min(k)) + Min(k);
%     g = G_hat(:,k);
%     GG = reshape(g,256,256);
%     G(:,:,k) = GG;
%     G(:,:,k) = G(:,:,k).*(Max(k)-Min(k)) + Min(k);
%     h = H_hat(:,k);
%     HH = reshape(h,256,256);
%     H(:,:,k) = HH;
%     H(:,:,k) = H(:,:,k).*(Max(k)-Min(k)) + Min(k);
%     x = X_hat(:,k);
%     XX = reshape(x,256,256);
%     X(:,:,k) = XX;
%     X(:,:,k) = X(:,:,k).*(Max(k)-Min(k)) + Min(k);
%     y = Y_hat(:,k);
%     YY = reshape(y,256,256);
%     Y(:,:,k) = YY;
%     Y(:,:,k) = Y(:,:,k).*(Max(k)-Min(k)) + Min(k);
%     z = Z_hat(:,k);
%     ZZ = reshape(z,256,256);
%     Z(:,:,k) = ZZ;
%     Z(:,:,k) = Z(:,:,k).*(Max(k)-Min(k)) + Min(k);
    
    Hyp_noise(:,:,k) = Hyp_noise(:,:,k).* (Max(k)-Min(k)) + Min(k);
end


%% image inpainting
% initialize random seed
% rng(2);
rng('default');

bandth1 = 48;
bandth2 = 92;

% u = im2double(imread([in 'discharge_print_128x128.png'])); 
u1 = A(:, :, bandth1);
u2 = F(:, :, bandth1);
u3 = A(:, :, bandth2);
u4 = F(:, :, bandth2);
[M,N,~] = size(u1);
% figure
% imshow(u)
% title('Original')

% Mask
% mask = im2double(imread([in 'maskthin.png']));
indm = mask>0;
% figure; imshow(mask); title('Mask')

% Conditioning points
w = 3; % thickness of the conditioning border
indc = get_conditioning_points(mask,w);
% figure; imshow(double(indc)); title('Conditioning points')

% Estimate Gaussian model on the whole image
[t1,m1] = estimate_adsn_model(u1,M,N);
[t2,m2] = estimate_adsn_model(u2,M,N);
[t3,m3] = estimate_adsn_model(u3,M,N);
[t4,m4] = estimate_adsn_model(u4,M,N);

% Inpainting
[v1,ke1,innov1] = kriging_inpainting(u1.*(1-indm),m1,t1,indm,indc);
[v2,ke2,innov2] = kriging_inpainting(u2.*(1-indm),m2,t2,indm,indc);
[v3,ke3,innov3] = kriging_inpainting(u3.*(1-indm),m3,t3,indm,indc);
[v4,ke4,innov4] = kriging_inpainting(u4.*(1-indm),m4,t4,indm,indc);

%% =====================================================================================

% MPSNR
PSNR_1 = zeros(1,band);
PSNR_2 = zeros(1,band);
PSNR_3 = zeros(1,band);
PSNR_4 = zeros(1,band);
PSNR_5 = zeros(1,band);
PSNR_6 = zeros(1,band);
PSNR_7 = zeros(1,band);
PSNR_8 = zeros(1,band);
PSNR_9 = zeros(1,band);
PSNR_10 = zeros(1,band);
PSNR_11 = zeros(1,band);
PSNR_12 = zeros(1,band);

for t = 1:band
%      [PSNR_1(t),MSE_1,MAXERR_1,L2RAT_1] = measerr(A(:,:,t),Hyp_actual(:,:,t));
%      [PSNR_2(t),MSE_2,MAXERR_2,L2RAT_2] = measerr(B(:,:,t),Hyp_actual(:,:,t));
%      [PSNR_3(t),MSE_3,MAXERR_3,L2RAT_3] = measerr(C(:,:,t),Hyp_actual(:,:,t));
%      [PSNR_4(t),MSE_4,MAXERR_4,L2RAT_4] = measerr(D(:,:,t),Hyp_actual(:,:,t));
%      [PSNR_5(t),MSE_5,MAXERR_5,L2RAT_5] = measerr(E(:,:,t),Hyp_actual(:,:,t));    
%      [PSNR_6(t),MSE_6,MAXERR_6,L2RAT_6] = measerr(F(:,:,t),Hyp_actual(:,:,t)); 
%      [PSNR_7(t),MSE_7,MAXERR_7,L2RAT_7] = measerr(G(:,:,t),Hyp_actual(:,:,t));
%      [PSNR_8(t),MSE_8,MAXERR_8,L2RAT_8] = measerr(H(:,:,t),Hyp_actual(:,:,t));
%      [PSNR_9(t),MSE_9,MAXERR_9,L2RAT_9] = measerr(X(:,:,t),Hyp_actual(:,:,t));
%      [PSNR_10(t),MSE_10,MAXERR_10,L2RAT_10] = measerr(Y(:,:,t),Hyp_actual(:,:,t));
%      [PSNR_11(t),MSE_11,MAXERR_11,L2RAT_11] = measerr(Z(:,:,t),Hyp_actual(:,:,t));
%      [PSNR_12(t),MSE_12,MAXERR_12,L2RAT_12] = measerr(Hyp_noise(:,:,t),Hyp_actual(:,:,t));


     PSNR_1(t) = psnr(A(:,:,t),Hyp_actual(:,:,t));
%      PSNR_2(t) = psnr(B(:,:,t),Hyp_actual(:,:,t));
%      PSNR_3(t) = psnr(C(:,:,t),Hyp_actual(:,:,t));
%      PSNR_4(t) = psnr(D(:,:,t),Hyp_actual(:,:,t));
%      PSNR_5(t) = psnr(E(:,:,t),Hyp_actual(:,:,t));    
     PSNR_6(t) = psnr(F(:,:,t),Hyp_actual(:,:,t)); 
%      PSNR_7(t) = psnr(G(:,:,t),Hyp_actual(:,:,t));
%      PSNR_8(t) = psnr(H(:,:,t),Hyp_actual(:,:,t));
%      PSNR_9(t) = psnr(X(:,:,t),Hyp_actual(:,:,t));
%      PSNR_10(t) = psnr(Y(:,:,t),Hyp_actual(:,:,t));
%      PSNR_11(t) = psnr(Z(:,:,t),Hyp_actual(:,:,t));
     PSNR_12(t) = psnr(Hyp_noise(:,:,t),Hyp_actual(:,:,t));
end
% psnr_1 = psnr(v1,Hyp_actual(:,:,30));
% psnr_2 = psnr(v2,Hyp_actual(:,:,30));
% psnr_3 = psnr(v3,Hyp_actual(:,:,125));
% psnr_4 = psnr(v4,Hyp_actual(:,:,125));

MPSNR_1 = mean(PSNR_1);
MPSNR_2 = mean(PSNR_2);
MPSNR_3 = mean(PSNR_3);
MPSNR_4 = mean(PSNR_4);
MPSNR_5 = mean(PSNR_5);
MPSNR_6 = mean(PSNR_6);
MPSNR_7 = mean(PSNR_7);
MPSNR_8 = mean(PSNR_8);
MPSNR_9 = mean(PSNR_9);
MPSNR_10 = mean(PSNR_10);
MPSNR_11 = mean(PSNR_11);
MPSNR_12 = mean(PSNR_12);

% MSSIM
Hyp_actual = double(Hyp_actual);
SSIM_1 = zeros(1,band);
SSIM_2 = zeros(1,band);
SSIM_3 = zeros(1,band);
SSIM_4 = zeros(1,band);
SSIM_5 = zeros(1,band);
SSIM_6 = zeros(1,band);
SSIM_7 = zeros(1,band);
SSIM_8 = zeros(1,band);
SSIM_9 = zeros(1,band);
SSIM_10 = zeros(1,band);
SSIM_11 = zeros(1,band);
SSIM_12 = zeros(1,band);

for p = 1:band
% SSIM_1(p) = SSIM(A(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_2(p) = SSIM(B(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_3(p) = SSIM(C(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_4(p) = SSIM(D(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_5(p) = SSIM(E(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_6(p) = SSIM(F(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_7(p) = SSIM(G(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_8(p) = SSIM(H(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_9(p) = SSIM(X(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_10(p) = SSIM(Y(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_11(p) = SSIM(Z(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),Max(p)-Min(p));
% SSIM_12(p) = SSIM(Hyp_noise(:,:,p),Hyp_actual(:,:,p),[0.01 0.03],fspecial('gaussian', 11, 1.5),1);


SSIM_1(p) = ssim(A(:,:,p),Hyp_actual(:,:,p));
% SSIM_2(p) = ssim(B(:,:,p),Hyp_actual(:,:,p));
% SSIM_3(p) = ssim(C(:,:,p),Hyp_actual(:,:,p));
% SSIM_4(p) = ssim(D(:,:,p),Hyp_actual(:,:,p));
% SSIM_5(p) = ssim(E(:,:,p),Hyp_actual(:,:,p));
SSIM_6(p) = ssim(F(:,:,p),Hyp_actual(:,:,p));
% SSIM_7(p) = ssim(G(:,:,p),Hyp_actual(:,:,p));
% SSIM_8(p) = ssim(H(:,:,p),Hyp_actual(:,:,p));
% SSIM_9(p) = ssim(X(:,:,p),Hyp_actual(:,:,p));
% SSIM_10(p) = ssim(Y(:,:,p),Hyp_actual(:,:,p));
% SSIM_11(p) = ssim(Z(:,:,p),Hyp_actual(:,:,p));
SSIM_12(p) = ssim(Hyp_noise(:,:,p),Hyp_actual(:,:,p));
end
% ssim_1 = ssim(v1,Hyp_actual(:,:,30));
% ssim_2 = ssim(v2,Hyp_actual(:,:,30));
% ssim_3 = ssim(v3,Hyp_actual(:,:,125));
% ssim_4 = ssim(v4,Hyp_actual(:,:,125));

MSSIM_1 = mean(SSIM_1);
MSSIM_2 = mean(SSIM_2);
MSSIM_3 = mean(SSIM_3);
MSSIM_4 = mean(SSIM_4);
MSSIM_5 = mean(SSIM_5);
MSSIM_6 = mean(SSIM_6);
MSSIM_7 = mean(SSIM_7);
MSSIM_8 = mean(SSIM_8);
MSSIM_9 = mean(SSIM_9);
MSSIM_10 = mean(SSIM_10);
MSSIM_11 = mean(SSIM_11);
MSSIM_12 = mean(SSIM_12);


figure(1);
imshow(Hyp_actual(:,:,bandth2),[]);
% figure(2);
% imshow(Hyp_actual_normalized(:,:,bandth),[]);
figure(3);
imshow(Hyp_noise(:,:,bandth2),[]);
figure(4);
imshow(A(:,:,bandth2),[]);
% figure(5);
% imshow(B(:,:,band),[]);
% figure(6);
% imshow(C(:,:,band),[]);
% figure(7);
% imshow(D(:,:,band),[]);
% figure(8);
% imshow(E(:,:,band),[]);
figure(9);
imshow(F(:,:,bandth2),[]);
% figure(10);
% imshow(G(:,:,band),[]);
% figure(11);
% imshow(H(:,:,band),[]);
figure(12);
imshow(v1,[]);
figure(13);
imshow(v2,[]);
figure(14);
imshow(v3,[]);
figure(15);
imshow(v4,[]);


% savefile = 'file.mat';
% save(savefile, 'MPSNR_1', 'MPSNR_2', 'MPSNR_3', 'MPSNR_4', 'MPSNR_5', 'MPSNR_6', 'MPSNR_7', 'MPSNR_8', 'MPSNR_9', 'MSSIM_1', 'MSSIM_2', 'MSSIM_3', 'MSSIM_4', 'MSSIM_5', 'MSSIM_6', 'MSSIM_7', 'MSSIM_8', 'MSSIM_9')