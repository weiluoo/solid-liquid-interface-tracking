%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script allows to reproduce the experiments of the paper
%   "Gaussian Texture Inpainting" 
%   (Bruno Galerne, Arthur Leclaire)
%   Preprint MAP5 n°2016-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% Input images folder
in = 'input/';
% Results folder
res = 'results/';
mkdir(res)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG : Teaser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = im2double(imread([in 'discharge_print_128x128.png'])); 
[M,N,~] = size(u);
figure
imshow(u)
title('Original')

% Mask
indm = im2double(imread([in 'masksiims.png']));
% figure; imshow(indm); title('Mask')

% Conditioning points
w = 3; % thickness of the conditioning border
indc = get_conditioning_points(indm,w);
%  figure; imshow(double(indc)); title('conditioning points')

% Estimate Gaussian model on the masked image
[t,m] = estimate_adsn_model(u.*(1-indm),M,N,indm); % masked image

% % check the texture model
% z = adsn_periodic(t,repmat(m,[M N 1]));
% figure
% imshow(z)
% title('ADSN model')

rng(2);
clear options; 
options.verb = 1;
options.ep = 1e-15;
options.imax = 1000;
options.normal = 1;
[v,ke,innov] = gaussian_inpainting(u.*(1-indm),m,t,indm,indc,options);

figure
imshow(u.*(1-indm));
title('Masked texture')

figure
imshow(v)
title('Inpainted')

figure
imshow(ke);
title('Kriging component')

figure
imshow(innov);
title('Innovation component')
 
% % save the results
fol = [res 'teaser/']; mkdir(fol)
imwrite(u,[fol 'original.png']);
maskedu = u.*(1-indm); maskedu(indm>0.5)=1; % imshow(maskedu);
imwrite(maskedu,[fol 'masked.png'])
imwrite(double(indc),[fol 'cond.png'])
imwrite(v,[fol 'inpainted.png'])
imwrite(ke,[fol 'kriging.png'])
imwrite(innov,[fol 'innovation.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG: Estimation of an ADSN model on a masked exemplar
%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'0000_2';'0001_1';'0070_4'; ...
    'discharge_print_128x128'; ...
    'fabric_2001';'fabric_2002';'fabric_2003'; ...
    'BarkPine0009_11_thumbhuge_128'; ...
    'BrickGroutless0039_2_thumbhuge_128'; ...
    'BrickOldDirty0009_2_thumbhuge_128'; ...
    'ConcreteBare0244_2_thumbhuge_128'; ...
    'ConcreteBare0280_39_thumbhuge_128'; ...
    'ConcreteFloors0068_2_thumbhuge_128'; ...
    'fabric_2004'; ...
    'Farmland0008_19_thumbhuge_128'; ...
    'Leather0002_2_thumbhuge_128'; ...
    'MarbleGranite0007_5_thumbhuge_128'; ...
    'MetalBare0131_thumbhuge_128'; ...
    'paint_2001'; ...
    'paint_2005'; ...
	'RockSediment0006_7_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128';
    'sand_2001';
    'sand_2003';
    'tissu_patchwork_128x128';
    };

close all

rng(1)        
mask = 1;
switch mask
    case 1
        masktype = 'square';
        ex = 14; %4, 9, 12 
        name = names{ex};
        u = im2double(imread([in name '.png']));
        [M,N,C] = size(u);
        x1 = 16; y1 = 16; x2 = 113; y2 = 113;
        indm = zeros(M,N,C);
        indm(y1:y2,x1:x2,:) = 1;        
    case 2
        masktype = 'excurs';
        ex = 22; %19, 22
        name = names{ex};
        u = im2double(imread([in name '.png']));
        [M,N,~] = size(u);
        % % excursion set of a Gaussian texture 
        indm = randn(M,N);
        indm = gaussian_convol(indm,4);
        indm = indm>0.007;
        indm = repmat(indm,[1 1 3]);
    case 3
        masktype = 'bern';
        ex = 20; % 6, 7, 12, 15, 16, 19, 20
        name = names{ex};
        u = im2double(imread([in name '.png']));
        p = 0.75; % percentage of masked pixels
        indm = rand(M,N)<p;
        indm = repmat(indm,[1 1 3]);
    case 4
        masktype = 'thin';
        ex = 9; % 6, 7, 9, 10, 12, 15, 16, 19, 20
        name = names{ex};
        u = im2double(imread([in name '.png']));
        indm = im2double(imread([in 'maskthin.png']));
        indm = (indm>0.5);
    case 5
        masktype = 'mix';
        ex = 19; % 5, 6, 15, 19, 25
        name = names{ex};
        u = im2double(imread([in name '.png']));
        indm = im2double(imread([in 'maskmix2.png']));
        indm = (indm>0.5);
end

maskedu = u.*(1-indm); maskedu(indm>0.5)=1; imshow(maskedu);
subplot(1,3,1)
imshow(maskedu)
title('Masked exemplar')

rng(2)
[t0,m0] = estimate_adsn_model(u);
v0 = adsn(t0,repmat(m0,[M N 1]));
subplot(1,3,2)
imshow(v0)
title('Oracle ADSN')

rng(2)
[t,m] = estimate_adsn_model(maskedu,M,N,indm);
v = adsn(t,repmat(m,[M N 1]));
subplot(1,3,3)
imshow(v)
title('Estimated ADSN')

% save the results
fol = [res 'covestimation/']; mkdir(fol)
fol = [fol masktype '/']; mkdir(fol)
imwrite(v0,[fol 'oracle.png']);
imwrite(v,[fol 'adsn.png']);
imwrite(maskedu,[fol 'maskedu.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG: Inpainting Gaussian textures with oracle model (regular masks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'0000_2';'0001_1';'0070_4'; ...
    'discharge_print_128x128'; ...
    'fabric_2001';'fabric_2002';'fabric_2003'; ...
    'BarkPine0009_11_thumbhuge_128'; ...
    'BrickGroutless0039_2_thumbhuge_128'; ...
    'BrickOldDirty0009_2_thumbhuge_128'; ...
    'ConcreteBare0244_2_thumbhuge_128'; ...
    'ConcreteBare0280_39_thumbhuge_128'; ...
    'ConcreteFloors0068_2_thumbhuge_128'; ...
    'fabric_2004'; ...
    'Farmland0008_19_thumbhuge_128'; ...
    'Leather0002_2_thumbhuge_128'; ...
    'MarbleGranite0007_5_thumbhuge_128'; ...
    'MetalBare0131_thumbhuge_128'; ...
    'paint_2001'; ...
    'paint_2005'; ...
	'RockSediment0006_7_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128';
    'sand_2001';
    'sand_2003';
    'tissu_patchwork_128x128';
    };

per = 0; % use periodic model

for ex = 3 % [3 8 10 12 13 20 22 25]
    
name = names{ex};
u = im2double(imread([in name '.png']));

[M,N,C] = size(u);
[t,m] = estimate_adsn_model(u);
rng(1);
if per
    u = adsn_periodic(t,repmat(m,[M N 1]));
else
    u = adsn(t,repmat(m,[M N 1]));
end
figure
subplot(2,2,1)
imshow(u)
title('Gaussian texture')

% medium square
x1 = 26; y1 = 26; x2 = 103; y2 = 103;
% big square
% x1 = 16; y1 = 16; x2 = 113; y2 = 113;
indm = zeros(M,N,C);
indm(y1:y2,x1:x2,:) = 1;
maskedu = u.*(1-indm); maskedu(indm>0.5)=1; imshow(maskedu);

clear options;
options.per = per;
options.verb = 1;

% conditioning with all the known values
indc = 1-indm; 
rng(3)
v = gaussian_inpainting(u.*(1-indm),m,t,indm,indc,options);

subplot(2,2,3)
imshow(v)
title('inpainted with all known values')
drawnow

% conditioning with a border
width = 3;
indc = get_conditioning_points(indm,width);
rng(3)
vb = gaussian_inpainting(u.*(1-indm),m,t,indm,indc,options);

subplot(2,2,2)
imshow(vb)
title('inpainted with bordering values')
drawnow

% conditioning with all the known values (cgg not on normal equations)
options.normal = 0;
indc = 1-indm; 
rng(3)
vd = gaussian_inpainting(u.*(1-indm),m,t,indm,indc,options);

subplot(2,2,4)
imshow(vd)
title('inpainted with CGD on direct equations')
drawnow

% save the results
if per
    fol = [res 'gausstex_sqhole_per/']; mkdir(fol)
else
    fol = [res 'gausstex_sqhole_nonper/']; mkdir(fol)
end
imwrite(u,[fol name '_adsn.png']);
imwrite(maskedu,[fol name '_masked.png'])
imwrite(v,[fol name '_inpainted.png'])
imwrite(vb,[fol name '_inpainted_w' num2str(width) '.png'])
imwrite(vd,[fol name '_inpainted_direct.png'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG: Inpainting Gaussian textures with oracle model (irregular masks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'0000_2';'0001_1';'0070_4'; ...
    'discharge_print_128x128'; ...
    'fabric_2001';'fabric_2002';'fabric_2003'; ...
    'BarkPine0009_11_thumbhuge_128'; ...
    'BrickGroutless0039_2_thumbhuge_128'; ...
    'BrickOldDirty0009_2_thumbhuge_128'; ... % 10
    'ConcreteBare0244_2_thumbhuge_128'; ...
    'ConcreteBare0280_39_thumbhuge_128'; ...
    'ConcreteFloors0068_2_thumbhuge_128'; ...
    'fabric_2004'; ...
    'Farmland0008_19_thumbhuge_128'; ...
    'Leather0002_2_thumbhuge_128'; ...
    'MarbleGranite0007_5_thumbhuge_128'; ...
    'MetalBare0131_thumbhuge_128'; ...
    'paint_2001'; ...
    'paint_2005'; ...
	'RockSediment0006_7_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128';
    'sand_2001';
    'sand_2003';
    'tissu_patchwork_128x128';
    };

per = 0; % use periodic model or not

for ex = 10 % [3 8 10 12 13 20 22 25]
    
name = names{ex};
u = im2double(imread([in name '.png']));

[M,N,~] = size(u);
[t,m] = estimate_adsn_model(u);
rng(1);
if per
    u = adsn_periodic(t,repmat(m,[M N 1]));
else
    u = adsn(t,repmat(m,[M N 1]));
end
figure
subplot(1,3,1)
imshow(u)
title('Gaussian texture')
drawnow

% two choices of masks
if 1
    % excursion set of a Gaussian texture
    choice = 'excurs';
    indm = randn(M,N);
    indm = gaussian_convol(indm,2);
    indm = indm>0.01;
    indm = repmat(indm,[1 1 3]);
else
    % Bernoulli mask
    choice = 'bern';
    p = 0.8; % percentage of masked pixels
    indm = rand(M,N)<p;
    indm = repmat(indm,[1 1 3]);
end

subplot(1,3,2)
maskedu = u.*(1-indm); maskedu(indm>0.5)=1; imshow(maskedu);
title('Masked texture')

clear options;
options.per = per;
options.verb = 2;
options.normal = 1;

% conditioning with all the known values
indc = 1-indm; 
rng(3)
v = gaussian_inpainting(u.*(1-indm),m,t,indm,indc,options);

subplot(1,3,3)
imshow(v)
title('inpainted with all known values')
drawnow
pause(1)

% save the results
if per
    fol = [res 'gausstex_' choice 'hole_per/']; mkdir(fol)
else
    fol = [res 'gausstex_' choice 'hole_nonper/']; mkdir(fol)
end
imwrite(u,[fol name '_adsn.png']);
imwrite(maskedu,[fol name '_masked.png'])
imwrite(v,[fol name '_inpainted.png'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG: Extrapolation - Oracle model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warning: very long

names = {'IMG_7407';'IMG_9875';'wall1021';'wall1022';'wall1023'; ...
    'wall1033';'wall1043';'wood1010_a'};
ex = 6;  % included in the paper: 6, 8

name = names{ex};
u = im2double(imread([in name '.png']));
big = 1;
if big
    M = 427; N = 641; C = 3;
else
    M = 256; N = 384; C = 3;
end
u = u(1:M,1:N,:);
[t,m] = estimate_adsn_model(u);
rng(1)
u = adsn(t,repmat(m,[M N 1]));
imshow(u)

choice = 'cmla';
indc = im2double(imread([in 'mask' choice '_' num2str(N) 'x' num2str(M) '.png']));
indc = double(indc>0);
indm = 1-indc;
% figure
% imshow(indm);
% figure
% imshow(indc)

upretty = ones(M,N,C);
upretty(~indm) = u(~indm);
figure
imshow(upretty)
title('Input')
%

clear options
options.normal = 0;
options.imax = 1000;
options.ep = 1e-15;
options.verb = 1;
rng(3)
v = gaussian_inpainting(u.*(1-indm),m,t,indm,indc,options);

% % restrict to a conditioning border
width = 3; % border thickness
indc = get_conditioning_points(indm,width);
rng(3)
vb = gaussian_inpainting(u.*(1-indm),m,t,indm,indc,options);

% Comparison with synthesizing the texture on the mask
w = adsn(t,repmat(m,[M N 1]));
w(~indm) = u(~indm);

figure
imshow(upretty);
title('Masked texture')

figure
imshow(v);
title('Inpainted')

figure
imshow(vb);
title('Inpainted with a border conditioning')

figure
imshow(w);
title('Baseline comparison')

% save the results
if big
    fol = [res 'extrapolationbig_' choice '/']; mkdir(fol)
else
    fol = [res 'extrapolation_' choice '/']; mkdir(fol)
end
fol = [fol name '_'];
imwrite(upretty,[fol 'original.png'])
imwrite(v,[fol 'extrapolated.png'])
imwrite(vb,[fol 'extrapolated_w' num2str(width) '.png'])
imwrite(w,[fol 'baseline.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG: inpainting with a Gaussian model estimated outside the mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'0000_2';'0001_1';'0070_4'; ...
    'discharge_print_128x128'; ...
    'fabric_2001';'fabric_2002';'fabric_2003'; ...
    'BarkPine0009_11_thumbhuge_128'; ...
    'BrickGroutless0039_2_thumbhuge_128'; ...
    'BrickOldDirty0009_2_thumbhuge_128'; ... % 10
    'ConcreteBare0244_2_thumbhuge_128'; ...
    'ConcreteBare0280_39_thumbhuge_128'; ...
    'ConcreteFloors0068_2_thumbhuge_128'; ...
    'fabric_2004'; ...
    'Farmland0008_19_thumbhuge_128'; ...
    'Leather0002_2_thumbhuge_128'; ...
    'MarbleGranite0007_5_thumbhuge_128'; ...
    'MetalBare0131_thumbhuge_128'; ...
    'paint_2001'; ...
    'paint_2005'; ...
	'RockSediment0006_7_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128';
    'sand_2001';
    'sand_2003';
    'tissu_patchwork_128x128';
    };

for ex = 3:5
name = names{ex};
u = im2double(imread([in name '.png']));

[M,N,C] = size(u);
%figure
imshow(u)
title('Original')

% Creation of the mask
indmask = 2; % 1 rectangle, 2 disc
switch indmask
    case 1
        x1 = 20; y1 = 32;
        x2 = 60; y2 = 96;
        mask = zeros(M,N);
        mask(y1:y2,x1:x2) = 1;
        % estimate the Gaussian model outside the mask
        [t,m] = estimate_adsn_model(u(:,x2+1:N,:),M,N);
    case 2
        x = 35; r = 30; 
        mask = disc(M,N,64,x,r,1);
        % estimate the Gaussian model outside the mask
        [t,m] = estimate_adsn_model(u(:,x+r+1:N,:),M,N);
        uw = draw_rectangle(u,x+r+1,N,1,M,1);
end
indm = repmat(mask>0,[1 1 C]);

% Conditioning points
indc = get_conditioning_points(indm,3);
% figure; imshow(double(indc)); title('Conditioning points')

% Automatic estimator of ADSN model
[tbis,mbis] = estimate_adsn_model(u,M,N,indm);

rdseed = 2;
rng(rdseed)
v = gaussian_inpainting(u.*(1-indm),m,t,indm,indc);
rng(rdseed)
vbis = gaussian_inpainting(u.*(1-indm),mbis,tbis,indm,indc);

%figure
imshow(u.*(1-indm))
title('Masked texture')

%figure
imshow(v)
title('Inpainted')

%figure
imshow(vbis)
title('Inpainted bis')

% save the results
fol = [res 'examples/']; mkdir(fol);
imwrite(u,[fol name '_original.png']);
maskedu = uw.*(1-indm); maskedu(indm>0)=1; % imshow(maskedu);
imwrite(maskedu,[fol name '_masked.png']);
imwrite(v,[fol name '_inpainted.png'])
imwrite(vbis,[fol name '_inpaintedbis.png'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG: Examples of textural inpainting
%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'0000_2';'0001_1';'0070_4'; ...
    'discharge_print_128x128'; ...
    'fabric_2001';'fabric_2002';'fabric_2003'; ...
    'BarkPine0009_11_thumbhuge_128'; ...
    'BrickGroutless0039_2_thumbhuge_128'; ...
    'BrickOldDirty0009_2_thumbhuge_128'; ...
    'ConcreteBare0244_2_thumbhuge_128'; ...
    'ConcreteBare0280_39_thumbhuge_128'; ...
    'ConcreteFloors0068_2_thumbhuge_128'; ...
    'fabric_2004'; ...
    'Farmland0008_19_thumbhuge_128'; ...
    'Leather0002_2_thumbhuge_128'; ...
    'MarbleGranite0007_5_thumbhuge_128'; ...
    'MetalBare0131_thumbhuge_128'; ...
    'paint_2001'; ...
    'paint_2005'; ...
	'RockSediment0006_7_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128'; ...
    'sand_2001'; ...
    'sand_2003'; ...
    'tissu_patchwork_128x128'; ...
    };

% macrotextures
namesmacro = {'BrickOldDirty0009_2_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128'; ...
    'blue_fabric2_crop'; 'fabric1010_crop'; ...
    'FloorsRounded0041_2_thumbhuge_128'; ...
    'raisins_verts_128'; 'rera1008_unzoom4_crop';
    'wall1022';'wall1023';'wood1010_a';
    'Farmland0008_19_thumbhuge_128' ...
    };

% %% see all exemplars
% 
% for ex=1:length(namesmacro)
%     name = namesmacro{ex};
%     u = im2double(imread([in name '.png']));
%     imshow(u)
%     title([name ' , ex no ' num2str(ex)]) 
%     drawnow
%     pause(1)
% end

%

close all

rng(1)        
for example = 1:10;
switch example
    case 1
        masktype = 'square';
        ex = 12; %4, 9, 12, 14 
        name = names{ex};
        u = im2double(imread([in name '.png']));
        x1 = 16; y1 = 16; x2 = 113; y2 = 113;
        indm = zeros(size(u));
        indm(y1:y2,x1:x2,:) = 1;        
    case 2
        masktype = 'excurs';
        ex = 22; %19, 22
        name = names{ex};
        u = im2double(imread([in name '.png']));
        [M,N,C] = size(u);
        % % excursion set of a Gaussian texture 
        indm = randn(M,N);
        indm = gaussian_convol(indm,4);
        indm = indm>0.007;
        indm = repmat(indm,[1 1 C]);
    case 3
        masktype = 'bern';
        ex = 5; % 6, 7, 12, 15, 16, 19, 20
        name = names{ex};
        u = im2double(imread([in name '.png']));
        [M,N,C] = size(u);
        p = 0.75; % percentage of masked pixels
        indm = rand(M,N)<p;
        indm = repmat(indm,[1 1 C]);
    case 4
        masktype = 'thin';
        ex = 20; % 6+, 7, 9+, 12, 15, 16, 19, 20
        name = names{ex};
        u = im2double(imread([in name '.png']));
        indm = im2double(imread([in 'maskthin.png']));
        indm = (indm>0.5);
    case 5
        masktype = 'mix';
        ex = 6; % 5, 6++, 15, 19+, 25
        name = names{ex};
        u = im2double(imread([in name '.png']));
        indm = im2double(imread([in 'maskmix2.png']));
        indm = (indm>0.5);
    case 6
        masktype = 'square';
        ex = 1; 
        name = namesmacro{ex};
        u = im2double(imread([in name '.png']));
        x1 = 16; y1 = 16; x2 = 113; y2 = 113;
        indm = zeros(size(u));
        indm(y1:y2,x1:x2,:) = 1;        
    case 7
        masktype = 'excurs';
        ex = 3;
        name = namesmacro{ex};
        u = im2double(imread([in name '.png']));
        u = u(1:128,1:128,:);
        [M,N,C] = size(u);
        % % excursion set of a Gaussian texture 
        indm = randn(M,N);
        indm = gaussian_convol(indm,4);
        indm = indm>0.007;
        indm = repmat(indm,[1 1 C]);
    case 8
        masktype = 'bern';
        ex = 8;
        name = namesmacro{ex};
        u = im2double(imread([in name '.png']));
        u = u(1:128,1:128,:);
        [M,N,C] = size(u);
        p = 0.75; % percentage of masked pixels
        indm = rand(M,N)<p;
        indm = repmat(indm,[1 1 C]);
    case 9
        masktype = 'thin';
        ex = 10; 
        name = namesmacro{ex};
        u = im2double(imread([in name '.png']));
        u = u(1:128,1:128,:);
        indm = im2double(imread([in 'maskthin.png']));
        indm = (indm>0.5);
    case 10
        masktype = 'mix';
        ex = 4;
        name = namesmacro{ex};
        u = im2double(imread([in name '.png']));
        u = u(1:128,1:128,:);
        indm = im2double(imread([in 'maskmix2.png']));
        indm = (indm>0.5);
end

[M,N,~] = size(u);
        
maskedu = u.*(1-indm); maskedu(indm>0.5)=1; imshow(maskedu);

% oracle adsn
rng(2)
[t0,m0] = estimate_adsn_model(u);
v0 = adsn(t0,repmat(m0,[M N 1]));

% inpainting
indc = get_conditioning_points(indm,3);
rng(2)
[t,m] = estimate_adsn_model(u,M,N,indm);
w = gaussian_inpainting(maskedu,m,t,indm,indc);

figure
subplot(1,3,1)
imshow(maskedu)
title('Masked exemplar')
subplot(1,3,2)
imshow(v0)
title('Oracle ADSN')
subplot(1,3,3)
imshow(w)
title('Inpainted')
drawnow

% save
fol = [res 'examples2/']; mkdir(fol)
fol = [fol name '/']; mkdir(fol)
imwrite(v0,[fol 'oracle.png']);
imwrite(maskedu,[fol 'maskedu.png'])
imwrite(w,[fol 'inpainted.png'])
imwrite(double(indm),[fol 'mask.png'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG: Inpainting textural parts of an image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = im2double(imread([in 'P1010203_reduced.png']));
[M,N,C] = size(u);

imshow(u)

% Mask creation
indm = zeros(M,N);
indm(124:155,98:160) = 1;
indm(205:241,121:202) = 1;
indm(323:361,167:246) = 1;
indm(472:515,217:311) = 1;
indm(92:136,271:307) = 1;
indm(206:253,392:430) = 1;
indm(123:167,440:476) = 1;
indm(116:150,512:586) = 1;
indm(194:230,572:636) = 1;
indm(308:353,663:746) = 1;

indm = repmat(indm,[1 1 C]);
indc = get_conditioning_points(indm,3);
% figure; imshow(cond); title('Conditioning points');

maskedu = u;
maskedu(indm>0) =1;

% learn the Gaussian model
x1 = 300; x2 = 640;
y1 = 270; y2 = 456;
[t,m] = estimate_adsn_model(u(y1:y2,x1:x2,:),M,N);
% Visualize u with the used window
uw = draw_rectangle(maskedu,x1,x2,y1,y2,1);

% figure
% imshow(uw.*(1-maskrgb))
% title('original')
%
options.verb = 1;
options.imax = 100;
v = gaussian_inpainting(maskedu,m,t,indm,indc,options);

figure
imshow(v)
title('inpainted')

% save the results
uc = u(180:540,80:500,:);
vc = v(180:540,80:500,:);
fol = [res 'nancy/']; mkdir(fol)
imwrite(u,[fol 'original.png'])
imwrite(v,[fol 'inpainted.png'])
imwrite(uw,[fol 'maskedu.png'])
imwrite(uc,[fol 'original_crop.png'])
imwrite(vc,[fol 'inpainted_crop.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG: Impact of the conditioning border - quantitative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

names = {'0000_2';'0001_1';'0070_4'; ...
    'discharge_print_128x128'; ...
    'fabric_2001';'fabric_2002';'fabric_2003'; ...
    'BarkPine0009_11_thumbhuge_128'; ...
    'BrickGroutless0039_2_thumbhuge_128'; ...
    'BrickOldDirty0009_2_thumbhuge_128'; ...
    'ConcreteBare0244_2_thumbhuge_128'; ...
    'ConcreteBare0280_39_thumbhuge_128'; ...
    'ConcreteFloors0068_2_thumbhuge_128'; ...
    'fabric_2004'; ...
    'Farmland0008_19_thumbhuge_128'; ...
    'Leather0002_2_thumbhuge_128'; ...
    'MarbleGranite0007_5_thumbhuge_128'; ...
    'MetalBare0131_thumbhuge_128'; ...
    'paint_2001'; ...
    'paint_2005'; ...
	'RockSediment0006_7_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128';
    'sand_2001';
    'sand_2003';
    'tissu_patchwork_128x128';
    };

ex = 4;
name = names{ex};
u = rgb2gray(im2double(imread([in name '.png'])));
su = std(u(:));

[M,N,C] = size(u);
figure
imshow(u)
title('Original')
%
% WARNING: VERY LONG SIMULATION. Results are saved.
%   Change 0 in 1 if you want to launch this experiment.
if 0
[t,m] = estimate_adsn_model(u,M,N);

ft = fft2(t);
max(abs(ft(:)).^2)

% Mask
wm = 15; % half length of the mask
indm = zeros(M,N,C);
indm(63+(-wm:wm),63+(-wm:wm),:) = 1;
mn = sum(sum(indm(:,:,1))); % number of masked pixels
figure; imshow(indm); title('Mask')

% Reference Gaussian model
% wref = 20;
% indcref = get_conditioning_points(indm,wref);
indcref = 1-indm;
cn = sum(sum(indcref(:,:,1)));
[kcompref,condcovref,kcoeffref,kcompvisuref,~,GammaMixref,GammaCondref,GammaMask] = get_conditional_model(u.*(1-indm),m,t,indm,indcref);

% sqrt(trace(condcovref)/mn)

% check the kriging coeff
disp(['Linf error in the kriging system ' num2str(max(max(kcoeffref*GammaCondref - GammaMixref')))])
% check the global covariance with approximate kriging coeff:
A = kcoeffref*GammaMixref;
GammaMaskApprox = GammaMask - A - A' + 2*kcoeffref*GammaCondref*kcoeffref';% 
disp(['Linf distance between covariance matrices = ' num2str(max(max(GammaMaskApprox(:)-GammaMask(:))))]);
% % OT distance between covariance matrices
% % d = compute_dist_covmatrix(GammaMask,GammaMaskApprox);
% % disp(['Covariance on mask error = ' num2str(d)])

% Parameter ep for studying the
% impact of an error of the pseudo-inverse computation
% Compute the required residual error to get an output
%   error which is < ep*std(u)
ep = 1e-2;

wmax = 20;
dmean = zeros(wmax,1);
dcov = zeros(wmax,1);
dcov1 = zeros(wmax,1);
dmodel = zeros(wmax,1);
dmodel1 = zeros(wmax,1);
err2tab = zeros(wmax,1);
errinftab = zeros(wmax,1);
condtab = zeros(wmax,1);
for w=1:wmax
    disp(w)
    indc = get_conditioning_points(indm,w);
    [kcomp,condcov,kcoeff,kcompvisu,~,GammaMix,GammaCond] = get_conditional_model(u.*(1-indm),m,t,indm,indc);
    dmean(w) = norm(kcomp(:)-kcompref(:))/sqrt(C*mn);
    dcov(w) = compute_dist_covmatrix(condcov,condcovref); % Frobenius distance
    dcov1(w) = compute_dist_covmatrix(condcov,condcovref,1); % opt transport distance
    dmodel(w) = sqrt(dmean(w)^2 + dcov(w)^2);
    dmodel1(w) = sqrt(dmean(w)^2 + dcov1(w)^2);
    err2tab(w) = ep*su*sqrt(C*mn)/norm(GammaMix);
    errinftab(w) = ep*su/norm(GammaMix,Inf);
    condtab(w) = cond(GammaCond);
end

end
%
load([in 'quantitative.mat'])

set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesLineWidth',1)
set(0,'defaultAxesFontSize',14)

figall = figure;
plot(1:wmax,dmodel/su,'black',1:wmax,dmean/su,'r--',1:wmax,dcov/su,'b:')
xlabel('w')
ylabel('\times \sigma_u')
legend({'Gaussian Model','Mean value','Covariance'})
title('Distance to reference Gaussian model')

figmean = figure;
plot(dmean/su)
title('Distance between means')
xlabel('w')
ylabel('\times \sigma_u')
figcov = figure;
plot(1:wmax,dcov/su,'r',1:wmax,dcov1/su,'b')
title('Distance between covariances')
xlabel('w')
ylabel('\times \sigma_u')
legend({'Fro','Opt'})
figmodel = figure;
plot(1:wmax,dmodel/su,'r',1:wmax,dmodel1/su,'b')
title('Distance between models')
xlabel('w')
ylabel('\times \sigma_u')
legend({'Fro','Opt'})
figcond = figure;
plot(condtab)
title('Condition number')
xlabel('w')
% ylabel('\times \sigma_u')

% figure
% plot(err2tab)
% title('L2 Stopping criterion for pseudo-inversion')
% figure
% plot(errinftab)
% title('Linf Stopping criterion for pseudo-inversion')

% figure
% imshow(kcompvisu)
% title('kriging')

% %% with a random model on the pseudo-inverse error
% % PROBLEM: the error matrix must be a SPD matrix which is not guaranteed
% % with the following simulation
% 
% err = 1e-3;
% % error on the pseudo-inverse matrix
% errmat = err*(2*rand(C*cn)-1);
% % error variance 
% errvar = GammaMix'*errmat*GammaMix;
% % output marginal error
% errout = sqrt(trace(errvar)/(C*mn));
% 
% disp(['std of marginal error / std(u) = ' num2str(errout/su)])

% illustrate with masked original
masku = u; masku(indm>0) = 1;
maskucond = repmat(masku,[1 1 3]);
indgb = zeros(M,N,C); indgb(:,:,2) = indcref; indgb(:,:,3) = indcref;
maskucond(indgb>0) = 0;

% save the results (need to download export_fig)
fol = [res 'condsize/' name '/']; mkdir(fol)
% export_fig(figall,[fol 'dall.pdf'],'-pdf','-transparent')
% export_fig(figmean, [fol 'dmean.pdf'],'-pdf','-transparent')
% export_fig(figcov, [fol 'dcov.pdf'],'-pdf','-transparent')
% export_fig(figmodel, [fol 'dmodel.pdf'],'-pdf','-transparent')
% export_fig(figcond,[fol 'cond.pdf'],'-pdf','-transparent')
imwrite(masku,[fol 'masku.png'])
imwrite(maskucond,[fol 'maskucond.png'])
save([fol 'save.mat'],'wmax','dmodel','dmodel1','su','dmean','dcov','dcov1','condtab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIG: Impact of the conditioning border - qualitative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

names = {'0000_2';'0001_1';'0070_4'; ...
    'discharge_print_128x128'; ...
    'fabric_2001';'fabric_2002';'fabric_2003'; ...
    'BarkPine0009_11_thumbhuge_128'; ...
    'BrickGroutless0039_2_thumbhuge_128'; ...
    'BrickOldDirty0009_2_thumbhuge_128'; ...
    'ConcreteBare0244_2_thumbhuge_128'; ...
    'ConcreteBare0280_39_thumbhuge_128'; ...
    'ConcreteFloors0068_2_thumbhuge_128'; ...
    'fabric_2004'; ...
    'Farmland0008_19_thumbhuge_128'; ...
    'Leather0002_2_thumbhuge_128'; ...
    'MarbleGranite0007_5_thumbhuge_128'; ...
    'MetalBare0131_thumbhuge_128'; ...
    'paint_2001'; ...
    'paint_2005'; ...
	'RockSediment0006_7_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128';
    'sand_2001';
    'sand_2003';
    'tissu_patchwork_128x128';
    };

ex = 4;
name = names{ex};
fol = [res 'condsizequal/' name '/']; mkdir(fol)

u = im2double(imread([in name '.png']));
su = std(u(:));

[M,N,C] = size(u);
figure
imshow(u)
title('Original')
%
[t,m] = estimate_adsn_model(u,M,N);

% Mask
wm = 30; % half length of the mask
indm = zeros(M,N,C);
indm(63+(-wm:wm),63+(-wm:wm),:) = 1;
mn = sum(sum(indm(:,:,1))); % number of masked pixels
figure; imshow(indm); title('Mask')

clear options
options.ep = 1e-6;
options.verb = 0;
options.imax = 1e3;
indcref = 1-indm;

rng(1)
[vref,keref,innovref] = gaussian_inpainting(u.*(1-indm),m,t,indm,indcref,options);

figure
subplot(2,4,1)
imshow(vref)
title('inpainted')
subplot(2,4,5)
imshow(keref)
title('krig comp')

imwrite(u.*(1-indm),[fol 'maskedu.png'])
imwrite(vref,[fol 'inpainted.png'])
imwrite(keref,[fol 'kcomp.png'])

k = 0;
for w=[1 3 5]
    k = k+1;
    disp(['w = ' num2str(w)])
    indc = get_conditioning_points(indm,w);
    rng(1)
    [v,ke,innov] = gaussian_inpainting(u.*(1-indm),m,t,indm,indc,options);
    disp(['Distance Mean values = ' num2str(norm(ke(:)-keref(:))/(su*sqrt(C*mn)))]);
    disp(['Distance Inp Results = ' num2str(norm(v(:)-vref(:))/(su*sqrt(C*mn)))]);

    subplot(2,4,1+k)
    imshow(v)
    title(['w=' num2str(w)])
    subplot(2,4,5+k)
    imshow(ke)
    title(['w=' num2str(w)])

    imwrite(v,[fol 'inpaintedw' num2str(w) '.png'])
    imwrite(ke,[fol 'kcompw' num2str(w) '.png'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig : Comparison with Criminisi and TV inpainting
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize random seed
rng(2)

u = im2double(imread([in 'criminisi_original.png']));
rescrim = im2double(imread([in 'criminisi_inpainted.png']));
u = u(110:250,164:340,:);
rescrim = rescrim(110:250,164:340,:);
mask = im2double(imread([in 'criminisi_mask4.png'])); mask = (mask(:,:,1)>0);

figure
imshow(u)
[M,N,C] = size(u);

indm = repmat(mask,[1 1 C]);
%figure; imshow(mask); title('Mask');

% get the conditioning points
indc = get_conditioning_points(mask,3);
% delete the part on the road
% on the crop
indc(110:end,:) = 0;
% figure; imshow(double(indc)); title('Conditioning points');

% estimate the Gaussian model
x1 = 70; x2 = N;
y1 = 10; y2 = 60;
[t,m] = estimate_adsn_model(u(y1:y2,x1:x2,:),M,N);

% Visualize u with the estimation window
uw = draw_rectangle(u,x1,x2,y1,y2,1);

% Inpainting
[v,kvisu,~] = gaussian_inpainting(u.*(1-indm),m,t,indm,indc);

figure
imshow(u)
figure
imshow(v)
title('inpainted')
figure
imshow(kvisu)
title('kriging component')

% save the results
fol = [res 'comp/']; mkdir(fol);
imwrite(u,[fol 'original.png'])
imwrite(uw,[fol 'original_window.png'])
imwrite(v,[fol 'inpainted.png'])
imwrite(rescrim,[fol 'criminisi.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig: Comparison with patch-based methods (wood texture)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = 'wood1008_sq_a_reduced';
u = im2double(imread([in name '.png']));
[M,N,~] = size(u);
meanu = mean(mean(u,2));
% figure
% imshow(u)
% title('Original')

% Creation of the mask
name = 'wood1008_sq_a_reduced_bis';
x1 = 40; y1 = 128;
x2 = 216; y2 = 230;
mask = zeros(M,N);
mask(y1:y2,x1:x2) = 1;
indm = repmat(mask,[1 1 3]);
[t,m] = estimate_adsn_model(u,M,N,indm);

% manual choice for estimating window
% xo1 = 1; yo1 = 1;
% xo2 = N; yo2 = y1-1;
% [t,m] = estimate_adsn_model(u(yo1:yo2,xo1:xo2,:),M,N);
% uw = draw_rectangle(u,xo1,xo2,yo1,yo2,2);
% figure; imshow(mask); title('Mask')
% imshow(uw.*(1-indm))

% Conditioning points
indc = get_conditioning_points(mask,3);
% figure; imshow(double(indc)); title('Conditioning points')

% Inpainting
rng(5)
options.verb = 1;
options.imax = 100;
[v,kc,innov] = gaussian_inpainting(u.*(1-indm),m,t,indm,indc,options);
v = v(1:M,1:N,:); kc = kc(1:M,1:N,:); innov = innov(1:M,1:N,:);

figure
imshow(u.*(1-indm));
title('Masked texture')

figure
imshow(v);
title('Inpainted')

figure
imshow(kc);
title('Kriging component')

fol = [res 'comp2/']; mkdir(fol);
fol = [fol name '/']; mkdir(fol);
imwrite(u,[fol 'original.png'])
maskedu = u.*(1-indm); maskedu(indm>0.5)=1; % imshow(maskedu);
imwrite(maskedu,[fol 'masked.png'])
imwrite(v,[fol 'inpainted.png'])
imwrite(kc,[fol 'kc.png'])
imwrite(innov,[fol 'innov.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IPOL PAPER Fig. CONVERGENCE OF THE CONJUGATE GRADIENT METHOD. SIMPLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('export_fig')

names = {'0000_2';'0001_1';'0070_4'; ...
    'discharge_print_128x128'; ...
    'fabric_2001';'fabric_2002';'fabric_2003'; ...
    'BarkPine0009_11_thumbhuge_128'; ...
    'BrickGroutless0039_2_thumbhuge_128'; ...
    'BrickOldDirty0009_2_thumbhuge_128'; ... % 10
    'ConcreteBare0244_2_thumbhuge_128'; ...
    'ConcreteBare0280_39_thumbhuge_128'; ...
    'ConcreteFloors0068_2_thumbhuge_128'; ...
    'fabric_2004'; ...
    'Farmland0008_19_thumbhuge_128'; ...
    'Leather0002_2_thumbhuge_128'; ...
    'MarbleGranite0007_5_thumbhuge_128'; ...
    'MetalBare0131_thumbhuge_128'; ...
    'paint_2001'; ...
    'paint_2005'; ...
	'RockSediment0006_7_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128';
    'sand_2001';
    'sand_2003';
    'tissu_patchwork_128x128';
    };

per = 0; % use periodic model or not

ex =  8; 
%[3 8 10 12 13 18 20 22 25]
% pretty: 1, 3, 5, 8, 10, 12, 13, 16, 18, 19, 20, 22, 25
% very pretty: 3, 8, 13, 18, 20, 22, 25
name = names{ex};
u = im2double(imread([in name '.png']));
% u = rgb2gray(u);

% crop
u = u(1:64,1:64,:);
[M,N,C] = size(u);

wm = 5; % half length of the mask
indm = zeros(size(u));
indm(32+(-wm:wm),32+(-wm:wm),:) = 1;

% [t,m] = estimate_adsn_model(u);
% u = adsn(t,repmat(m,[M N 1]));

% % get red model
% u(:,:,2:3) = 0;
% t(:,:,2:3) = 0;
% m(:,:,2:3) = 0;

% indc = 1-indm;
indc = get_conditioning_points(indm,3);
maskedu = u.*(1-indm); maskedu(indm>0.5)=1; 
figure
imshow(maskedu);
title('Masked texture')

[t,m] = estimate_adsn_model(maskedu,M,N,indm);

clear options;
options.per = per;
options.reg = 0; 
options.condnum = 1;

rdseed = 3;

rng(rdseed)
vref = kriging_inpainting(u,m,t,indm,indc,options);

options.vref = vref;
options.imax = 10000;
options.ep = 0;
options.verb = 1;
options.normal = 1;

rng(rdseed)
[v,distvref,distiter,resl2norm,reslinfnorm] = gaussian_inpainting_convspeed(u,m,t,indm,indc,options);

% with regularization
options.reg = 0.01; % add white noise of std reg in the model
rng(rdseed)
vrefreg = kriging_inpainting(u,m,t,indm,indc,options);
options.vref = vrefreg;
rng(rdseed)
[vreg,distvrefreg,distiterreg,resl2normreg,reslinfnormreg] = gaussian_inpainting_convspeed(u,m,t,indm,indc,options);

figure
imshow(v)
title('inpainted fast')
drawnow
figure
imshow(vref)    
title('inpainted ref')
drawnow
figure
imshow(vreg)
title('reg inpainted fast')
drawnow
figure
imshow(vrefreg)    
title('reg inpainted ref')
drawnow

set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesLineWidth',1)
set(0,'defaultAxesFontSize',12)

mn = sum(indm(:));

figdistvref = figure;
semilogy(distvref/sqrt(C*mn))
title('Inpainting L² error err(k)')
xlabel('Number of Iterations k')

figres = figure;
semilogy(resl2norm)
title('Residual L² norm resn(k)')
xlabel('Number of Iterations k')

figdistvrefreg = figure;
semilogy(distvrefreg/sqrt(C*mn))
title('Inpainting L² error err(k)')
xlabel('Number of Iterations k')

figresreg = figure;
semilogy(resl2normreg)
title('Residual L² norm resn(k)')
xlabel('Number of Iterations k')

% figure
% plot(distiter/sqrt(C*mn))
% title('Distance between iterates')
% figure
% semilogx(distiter/sqrt(C*mn))
% title('Distance between iterates log scale')

% save the results (need export_fig)
fol = [res 'convergence/']; mkdir(fol)
% export_fig(figdistvref, [fol 'distvref.pdf'],'-pdf','-transparent')
% export_fig(figres, [fol 'resl2norm.pdf'],'-pdf','-transparent')
% export_fig(figdistvrefreg, [fol 'distvrefreg.pdf'],'-pdf','-transparent')
% export_fig(figresreg, [fol 'resl2normreg.pdf'],'-pdf','-transparent')
imwrite(maskedu,[fol 'maskedu.png'])
imwrite(vref,[fol 'vref.png'])
imwrite(v,[fol 'v.png'])
imwrite(vrefreg,[fol 'vrefreg.png'])
imwrite(vreg,[fol 'vreg.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. IPOL PAPER CONVERGENCE OF THE CONJUGATE GRADIENT METHOD. 
% Difficult case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'0000_2';'0001_1';'0070_4'; ...
    'discharge_print_128x128'; ...
    'fabric_2001';'fabric_2002';'fabric_2003'; ...
    'BarkPine0009_11_thumbhuge_128'; ...
    'BrickGroutless0039_2_thumbhuge_128'; ...
    'BrickOldDirty0009_2_thumbhuge_128'; ... % 10
    'ConcreteBare0244_2_thumbhuge_128'; ...
    'ConcreteBare0280_39_thumbhuge_128'; ...
    'ConcreteFloors0068_2_thumbhuge_128'; ...
    'fabric_2004'; ...
    'Farmland0008_19_thumbhuge_128'; ...
    'Leather0002_2_thumbhuge_128'; ...
    'MarbleGranite0007_5_thumbhuge_128'; ...
    'MetalBare0131_thumbhuge_128'; ...
    'paint_2001'; ...
    'paint_2005'; ...
	'RockSediment0006_7_thumbhuge_128'; ...
    'RoofingPlates0017_1_thumbhuge_128';
    'sand_2001';
    'sand_2003';
    'tissu_patchwork_128x128';
    };

per = 0; % use periodic model or not

ex =  8; 
%[3 8 10 12 13 18 20 22 25]
% pretty: 1, 3, 5, 8, 10, 12, 13, 16, 18, 19, 20, 22, 25
% very pretty: 3, 8, 13, 18, 20, 22, 25
name = names{ex};
u = im2double(imread([in name '.png']));
% u = rgb2gray(u);

% crop
u = u(1:64,1:64,:);
[M,N,C] = size(u);

wm = 5; % half length of the mask
indm = zeros(size(u));
indm(32+(-wm:wm),32+(-wm:wm),:) = 1;

% [t,m] = estimate_adsn_model(u);
% u = adsn(t,repmat(m,[M N 1]));

% % get red model
% u(:,:,2:3) = 0;
% t(:,:,2:3) = 0;
% m(:,:,2:3) = 0;

indc = 1-indm;
% indc = get_conditioning_points(indm,3);
maskedu = u.*(1-indm); maskedu(indm>0.5)=1; 
figure
imshow(maskedu);
title('Masked texture')

[t,m] = estimate_adsn_model(maskedu,M,N,indm);

clear options;
options.per = per;
options.reg = 0; 
options.imax = 10000;
options.ep = 0;
options.verb = 1;
options.normal = 1;

rdseed = 3;

rng(rdseed)
[v,~,distiter,resl2norm,reslinfnorm] = gaussian_inpainting_convspeed(u,m,t,indm,indc,options);

% with regularization
options.reg = 0.01; % add white noise of std reg in the model
rng(rdseed)
[vreg,~,distiterreg,resl2normreg,reslinfnormreg] = gaussian_inpainting_convspeed(u,m,t,indm,indc,options);

figure
imshow(v)
title('inpainted fast')
drawnow
figure
imshow(vreg)
title('reg inpainted fast')
drawnow

set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesLineWidth',1)
set(0,'defaultAxesFontSize',12)

figres = figure;
semilogy(resl2norm)
title('Residual L² norm resn(k)')
xlabel('Number of Iterations k')

figresreg = figure;
semilogy(resl2normreg)
title('Residual L² norm resn(k)')
xlabel('Number of Iterations k')

% figure
% plot(distiter/sqrt(C*mn))
% title('Distance between iterates')
% figure
% semilogx(distiter/sqrt(C*mn))
% title('Distance between iterates log scale')

% save the results (need export_fig)
fol = [res 'convergence2/']; mkdir(fol)
% export_fig(figres, [fol 'resl2norm.pdf'],'-pdf','-transparent')
% export_fig(figresreg, [fol 'resl2normreg.pdf'],'-pdf','-transparent')
imwrite(maskedu,[fol 'maskedu.png'])
imwrite(v,[fol 'v.png'])
imwrite(vreg,[fol 'vreg.png'])



