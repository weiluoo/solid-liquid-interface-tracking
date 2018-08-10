function [v,kvisu,innov,precsys,csys] = kriging_inpainting(u,mv,t,indm,indc,options)
%KRIGING_INPAINTING  
%   Microtexture inpainting with Gaussian conditional simulation.
%
% [v,kvisu,innov,precsys,csys] = kriging_inpainting(u,mv,t,indm,indc,options)
%
%   Inpaint u in the masked region (given by indm) by drawing
%   a conditional sample of the ADSN model associated to the mean value mv
%   and the texton t, given values on a conditioning set (given by indc).
%
%   This software has been improved in gaussian_inpainting.m where
%    the linear systems are solved with a conjugate gradient descent.
%   
%   INPUT
%   u       Initial image
%   mv      Mean value of the estimated Gaussian model (ADSN)
%   t       Texton of the estimated Gaussian model (ADSN)
%   indm    Indicator function of masked region
%   indc    Indicator function of conditioning region
%   options Structure containing optional inputs:
%   - options.per  (default 0)
%       Use periodic ADSN if per==1, otherwise use non-periodic ADSN
%   - options.reg (default 0)
%       Add regularization in the linear system
%   - options.condnum (default 0)
%       Print condition number of the linear system
%
%   OUTPUT
%   v       Inpainted image
%   kvisu   (optional) Kriging component
%   innov   (optional) Innovation component
%   csys    (optional) Condition number of the linear system
%   precsys (optional) Precision on the solution
%
%   NB:
%   - The mean value is reimposed on kvisu and innov for visualization.
%   - indm can have dimensions MxNxC (provided that the
%      arrays indm(:,:,c) are the same for c=1:C). Idem for indc.
%
%   WARNING:
%   Computing the condition number of the system may take some time!
%   
%   See detailed description of the algorithm in the paper
%   "Microtexture Inpainting through Gaussian Conditional Simulation" 
%   (Bruno Galerne, Arthur Leclaire, Lionel Moisan)
%   Proceedings of the International Conference on Acoustics, Speech, and
%   Signal Processing, 2016.
%
%   If you use this code for further work, please be sure to make
%   reference to this paper.
%
%   Author: Arthur Leclaire
%   v 1.0 (03/2016) First public version
%   v 1.1 (04/2016) Corrected bug in covariance computation
%   v 1.2 (05/2016) Display residual norm
%   v 1.3 (05/2016) Added regularization parameter + condnum option

options.trash = 0;

if ~isfield(options,'per')
    per = 0;
else
    per = options.per;
end

if ~isfield(options,'reg')
    reg = 0;
else
    reg = options.reg;
end

if ~isfield(options,'condnum')
    condnum = (nargout>4);
else
    condnum = (options.condnum)||(nargout>4);
end

[M,N,C] = size(u);

% Handle the dimensions of indm and indc
indm = indm(:,:,1);
indc = indc(:,:,1);

% technical restriction
% (in order to simplify the covariance computation)
if (size(t,1)>M)||(size(t,2)>N)
    error('The texton must not be larger than the input');
end

if per==0
    M = 2*M; N = 2*N;
    u = zeropad(u,M,N);
    t = zeropad(t,M,N);
    indm = zeropad(indm,M,N);
    indc = zeropad(indc,M,N);
end

mp = find(indm>0);
cp = find(indc>0);
cn = length(cp);
[ca,cb] = ind2sub([M N],cp);

% Covariance of the fitted random field
covt = zeros(M,N,C,C);
ft = fft2(t);
for a=1:C
    for b=1:C
        covt(:,:,a,b) = real(ifft2(conj(ft(:,:,a)).*ft(:,:,b)));
    end
end

% Fill the matrix of the system
S = zeros(C*cn);
for j=1:cn
    for a=1:C
        for b=1:C
            shifta = cn*(a-1);
            shiftb = cn*(b-1);
            % covariance around (ca(j),cb(j))
            covtt = circshift(covt(:,:,a,b),[ca(j)-1 cb(j)-1]);
            % extract the part on the conditioning points
            S(shifta+j,shiftb+(1:cn)) = covtt(cp);
        end
    end
end
S = S + reg^2*eye(C*cn);

% initialize the output with the masked texture
v = u;

% draw an independent realization of the random field
z = adsn_periodic(t,zeros(M,N,C)) + reg*randn(M,N,C);
% NB: if per==1, t has already been extended by zero-padding
% so that the half-crop of z is indeed a sample of the
% non-periodic ADSN associated to t.

% if kvisu and innov are not required (i.e. nargout==1),
% then one only needs the kriging component of u-mv-z
%   (faster to compute)
if nargout>1
    
    % compute the RHS of the system
    T = zeros(C*cn,2);
    for c=1:C
        shiftc = cn*(c-1);
        uc = u(:,:,c); zc = z(:,:,c);
        T(shiftc+(1:cn),1) = uc(cp)-mv(:,:,c);
        T(shiftc+(1:cn),2) = zc(cp);
    end
    
    % system solve:
    tic;
    sol = S\T;
    timetosolve=toc;
    disp(['Time to solve the linear system of size ',num2str(C*cn),' x ' num2str(C*cn) ': ' num2str(timetosolve)]);
    r = T-S*sol;
    rp = S*r;
    rn2 = mean(r(:).^2); rn2 = sqrt(rn2);
    rpn2 = mean(rp(:).^2); rpn2 = sqrt(rpn2);    
    disp(['Residual L2 norm = ' num2str(rn2)])
    disp(['Normal Residual L2 norm = ' num2str(rpn2)])
    
    if nargout>3
        precsys = max(max(abs(S*sol-T)));
    end
    if condnum
        csys = cond(S);
        disp(['Condition number = ' num2str(csys)])
    end
    
    % compute the kriging components
    k = zeros(M,N,C); % kriging component
    zk = zeros(M,N,C); % kriging component of the ADSN
    for c=1:C
        shiftc = cn*(c-1);
        kc = zeros(M,N); zkc = zeros(M,N);
        kc(cp) = sol(shiftc+(1:cn),1);
        zkc(cp) = sol(shiftc+(1:cn),2);
        k(:,:,c) = kc;
        zk(:,:,c) = zkc;
    end
    k = conv_with_covariance(k,t,reg);
    zk = conv_with_covariance(zk,t,reg);
    
    kvisu = u;
    innov = repmat(mv,[M N 1]);
    for c=1:C
        kc = k(:,:,c); kvisuc = kvisu(:,:,c);
        zkc = zk(:,:,c); zc = z(:,:,c);
        vc = v(:,:,c);
        vc(mp) = mv(1,1,c) + kc(mp) + zc(mp) - zkc(mp);
        v(:,:,c) = vc;
        
        kvisuc(mp) = mv(1,1,c) + kc(mp);
        kvisu(:,:,c) = kvisuc;
        
        innovc = repmat(mv(1,1,c),[M N 1]);
        innovc(mp) = mv(1,1,c) + zc(mp) - zkc(mp);
        innov(:,:,c) = innovc;
    end
else
    
    % compute the RHS of the system
    T = zeros(C*cn,1);
    for c=1:C
        shiftc = cn*(c-1);
        uc = u(:,:,c); zc = z(:,:,c);
        T(shiftc+(1:cn),1) = uc(cp)-mv(:,:,c)-zc(cp);
    end
    
    % system solve:
    tic;
    sol = S\T;
    timetosolve=toc;
    disp(['Time to solve the linear system of size ',num2str(C*cn),' x ' num2str(C*cn) ': ' num2str(timetosolve)]);
    r = T-S*sol;
    rp = S*r;
    rn2 = mean(r(:).^2); rn2 = sqrt(rn2);
    rpn2 = mean(rp(:).^2); rpn2 = sqrt(rpn2);    
    disp(['Residual L2 norm = ' num2str(rn2)])
    disp(['Normal Residual L2 norm = ' num2str(rpn2)])
    
    if nargout>3
        precsys = max(max(abs(S*sol-T)));
    end
    if condnum
        csys = cond(S);
        disp(['Condition number = ' num2str(csys)])
    end
    
    % compute the kriging component u-mv-z
    k = zeros(M,N,C);
    for c=1:C
        shiftc = cn*(c-1);
        kc = zeros(M,N);
        kc(cp) = sol(shiftc+(1:cn),1);
        k(:,:,c) = kc;
    end
    k = conv_with_covariance(k,t,reg);
    
    for c=1:C
        kc = k(:,:,c);
        zc = z(:,:,c);
        vc = v(:,:,c);
        vc(mp) = mv(1,1,c) + kc(mp) + zc(mp);
        v(:,:,c) = vc;
    end
end

% If we used a non-periodic model, a half-crop is needed to come back
%   to the original image domain.
if per==0
    M = floor(M/2); N = floor(N/2);
    v = v(1:M,1:N,:);
    if nargout>1
        kvisu = kvisu(1:M,1:N,:);
        innov = innov(1:M,1:N,:);
    end
end
end

function v = conv_with_covariance(u,t,reg)
[M,N,C] = size(u);
tmp = zeros(M,N);
v = zeros(M,N,C);
ft = fft2(t);
for c=1:C
    tmp = tmp + conj(ft(:,:,c)).*fft2(u(:,:,c));
end
for c=1:C
    v(:,:,c) = real(ifft2(ft(:,:,c).*tmp));
end
v = v + reg^2*u;
end
