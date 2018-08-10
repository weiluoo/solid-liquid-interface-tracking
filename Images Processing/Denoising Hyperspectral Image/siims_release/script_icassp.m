
clear all
close all

% % Input images folder
% in = 'input/';
% % Results folder
% res = 'results/';
% mkdir(res)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig1 : Inpainting with oracle Gaussian model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize random seed
% rng(2);
rng('default');

% u = im2double(imread([in 'discharge_print_128x128.png'])); 
u = double(imread('C:/users/lh248/downloads/image2.jpg'));
u = u(:, :, 1);
[M,N,~] = size(u);
figure
imshow(u, [])
title('Original')

% Mask
% mask = im2double(imread([in 'maskthin.png']));
mask = zeros(M, N);
mask(100:101, 101:200) = 1;
indm = mask>0;
figure; imshow(mask); title('Mask')

% Conditioning points
w = 3; % thickness of the conditioning border
indc = get_conditioning_points(mask,w);
figure; imshow(double(indc)); title('Conditioning points')

% Estimate Gaussian model on the whole image
[t,m] = estimate_adsn_model(u,M,N);
% % check the texture model
% z = adsn_periodic(t,repmat(m,[M N 1]));
% figure
% imshow(z)
% title('ADSN model')

% Inpainting
[v,ke,innov] = kriging_inpainting(u.*(1-indm),m,t,indm,indc);

figure
imshow(u.*(1-mask), [])
title('Masked texture')

figure
imshow(v, [])
title('Inpainted')

figure
imshow(ke, [])
title('Kriging component')

figure
imshow(innov, [])
title('Innovation component')
 
% % save the results
% fol = [res 'validation/']; mkdir(fol)
% imwrite(u,[fol 'original.png']);
% maskedu = u.*(1-indm); maskedu(indm>0.5)=1; % imshow(maskedu);
% imwrite(maskedu,[fol 'masked.png'])
% imwrite(indc,[fol 'cond.png'])
% imwrite(v,[fol 'inpainted.png'])
% imwrite(ke,[fol 'kriging.png'])
% imwrite(innov,[fol 'innovation.png'])

toc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Fig2 : Inpainting with a Gaussian model leanrt outside the mask
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tic;
% 
% % initialize random seed
% % rng(2);
% rng('default');
% 
% names = {'0000_2';'0001_1';'0070_4'; ...
%     'discharge_print_128x128'; ...
%     'fabric_2001';'fabric_2002';'fabric_2003'};
% ex = 1;
% name = names{ex};
% u = im2double(imread([in name '.png']));
% 
% [M,N,C] = size(u);
% figure
% imshow(u)
% title('Original')
% 
% % Mask and estimation of the Gaussian model outside the mask
% indmask = 2; % 1 rectangle, 2 disc
% switch indmask
%     case 1
%         x1 = 20; y1 = 32;
%         x2 = 60; y2 = 96;
%         mask = zeros(M,N,C);
%         mask(y1:y2,x1:x2,:) = 1;
%         % Estimate the Gaussian model outside the mask
%         [t,m] = estimate_adsn_model(u(:,x2+1:N,:),M,N);
%         uw = draw_rectangle(u,x1,y1,x2,y2,1);
%     case 2
%         x = 35; r = 30; 
%         mask = repmat(disc(M,N,64,x,r,1),[1 1 C]);
%         % Estimate the Gaussian model outside the mask
%         [t,m] = estimate_adsn_model(u(:,x+r+1:N,:),M,N);
%         uw = draw_rectangle(u,x+r+1,N,1,M,1);
% end
% indm = mask>0;
% figure; imshow(mask); title('Mask')
% % % test the texture model
% % z = adsn_periodic(t,repmat(m,[M N 1]));
% % figure
% % imshow(z)
% % title('ADSN model')
% 
% % Conditioning points
% indc = get_conditioning_points(mask,3);
% figure; imshow(double(indc)); title('Conditioning points')
% 
% % Inpainting
% v = kriging_inpainting(u.*(1-indm),m,t,indm,indc);
% 
% figure
% imshow(u.*(1-mask))
% title('Masked texture')
% 
% figure
% imshow(v)
% title('Inpainted')
% 
% % % save the results
% % fol = [res 'examples/']; mkdir(fol);
% % imwrite(u,[fol name '_original.png']);
% % maskedu = uw.*(1-indm); maskedu(indm>0)=1; % imshow(maskedu);
% % imwrite(maskedu,[fol name '_masked.png']);
% % imwrite(v,[fol name '_inpainted.png'])
% 
% toc;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Fig3 : Comparison with Criminisi
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% % initialize random seed
% % rng(2)
% rng('default');
% 
% u = im2double(imread([in 'criminisi_original.png']));
% rescrim = im2double(imread([in 'criminisi_inpainted.png']));
% u = u(110:250,164:340,:);
% rescrim = rescrim(110:250,164:340,:);
% mask = im2double(imread([in 'criminisi_mask.png']));
% mask = repmat(mask(:,:,1),[1 1 3]);
% indm = mask>0;
% 
% figure
% imshow(u)
% [M,N,C] = size(u);
% 
% figure; imshow(mask); title('Mask');
% 
% % get the conditioning points
% indc = get_conditioning_points(mask,3);
% % delete the part on the road
% indc(110:end,:) = 0;
% figure; imshow(double(indc)); title('Conditioning points');
% 
% % Estimate the Gaussian model
% x1 = 70; x2 = N;
% y1 = 10; y2 = 80;
% [t,m] = estimate_adsn_model(u(y1:y2,x1:x2,:),M,N);
% uw = draw_rectangle(u,x1,x2,y1,y2,1);
% figure; imshow(uw); title('Red box for estimation of ADSN Model');
% 
% % %% test the texture model
% % z = adsn_periodic(t,repmat(m,[M N 1]));
% % figure
% % imshow(z)
% % title('ADSN model')
% 
% [v,kvisu,~] = kriging_inpainting(u.*(1-indm),m,t,indm,indc);
% 
% figure
% imshow(u)
% figure
% imshow(v)
% title('inpainted')
% figure
% imshow(kvisu)
% title('kriging component')
% 
% % fol = [res 'comp/']; mkdir(fol);
% % imwrite(u,[fol 'original.png'])
% % imwrite(uw,[fol 'original_window.png'])
% % imwrite(v,[fol 'inpainted.png'])
% % imwrite(rescrim,[fol 'criminisi.png'])
% 
% toc;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Fig4 : Comparison with Efros-Leung
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% 
% name = 'wood1008_sq_a_reduced';
% u = im2double(imread([in name '.png']));
% [M,N,C] = size(u);
% meanu = mean(mean(u,2));
% figure
% imshow(u)
% title('Original')
% 
% % Mask
% x1 = 40; y1 = 128;
% x2 = 216; y2 = 230;
% indm = zeros(M,N,C);
% indm(y1:y2,x1:x2,:) = 1;
% % Estimate the Gaussian model outside the mask
% xo1 = 1; yo1 = 1;
% xo2 = N; yo2 = y1-1;
% [t,m] = estimate_adsn_model(u(yo1:yo2,xo1:xo2,:),M,N);
% uw = draw_rectangle(u,xo1,xo2,yo1,yo2,2);
% figure
% imshow(uw.*(1-indm))
% title('Masked original with ADSN red box')
% % % test the texture model
% % z = adsn_periodic(t,repmat(m,[M N 1]));
% % figure
% % imshow(z)
% % title('ADSN model')
% 
% % Conditioning points
% indc = get_conditioning_points(indm,3);
% cp = find(indc>0);
% [ca,cb] = ind2sub([M N],cp);
% figure; imshow(double(indc)); title('Conditioning points')
% 
% 
% % rng(2)
% rng('default');
% [v,kc,innov] = kriging_inpainting(u.*(1-indm),m,t,indm,indc);
% 
% figure
% imshow(u.*(1-indm));
% title('Masked texture')
% 
% figure
% imshow(v);
% title('Inpainted')
% 
% figure
% imshow(kc);
% title('Kriging component')
% 
% % fol = [res 'comp2/']; mkdir(fol);
% % fol = [fol name '/']; mkdir(fol);
% % imwrite(u,[fol 'original.png'])
% % maskedu = uw.*(1-indm); maskedu(indm>0.5)=1; % imshow(maskedu);
% % imwrite(maskedu,[fol 'masked.png'])
% % imwrite(v,[fol 'inpainted.png'])
% % imwrite(kc,[fol 'kc.png'])
% % imwrite(innov,[fol 'innov.png'])
% toc;
