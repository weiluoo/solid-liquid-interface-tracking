%% Microtexture Inpainting by Gaussian Conditional Simulation - Tutorial
  
close all
clear all

in = 'input/';

rng(3); % initialize random seed

%% Load input image

name = 'wood1008_sq_a_reduced';
u = im2double(imread([in name '.png']));
[M,N,C] = size(u);
meanu = mean(mean(u,2));

figure
imshow(u)
title('Original')

%% Create a mask

x1 = 40; y1 = 128;
x2 = 216; y2 = 230;
indm = zeros(M,N,C);
indm(y1:y2,x1:x2,:) = 1;

figure
imshow(indm)
title('Mask')

%% Compute the Conditioning points on the mask border

indc = get_conditioning_points(indm,3);
figure
imshow(double(indc))
title('Conditioning points')

%% Estimate a Gaussian model outside the mask
xo1 = 1; yo1 = 1;
xo2 = N; yo2 = y1-1;
[t,m] = estimate_adsn_model(u(yo1:yo2,xo1:xo2,:),M,N);
uw = draw_rectangle(u,xo1,xo2,yo1,yo2,2);
figure
imshow(uw.*(1-indm))
title(sprintf('Masked original. ADSN estimated in the red box'))

%% Check the texture model

z = adsn_periodic(t,repmat(m,[M N 1]));

figure
imshow(z)
title('Realization of the ADSN model')

%% Inpaint by Gaussian conditional simulation

[v,kc,innov] = gaussian_inpainting(u.*(1-indm),m,t,indm,indc);

figure
imshow(u.*(1-indm));
title('Masked texture')
drawnow

figure
imshow(v);
title('Inpainted')
drawnow

figure
imshow(kc);
title('Kriging component')
drawnow

figure
imshow(innov);
title('Innovation component')
drawnow
